//! `aks`: run the Lenstra–AKS primality test on one number.
//!
//! A superset of the corrected C++ programs' command line:
//!
//! ```text
//! aks [--theorem | --phi | --lambda] [n]        # n from standard input if omitted
//!     [--mul kronecker|kronecker-bytes|schoolbook] [-j N] [--log FILE] [--quiet]
//! ```
//!
//! Exit status 0 for prime, 1 for composite (and for n < 2), 2 for bad usage,
//! as in the C++, so the same sweep scripts work with either.

use std::fs::OpenOptions;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use aks2::aks::{ARange, Multiplier, Options, Reason, Trace, aks};
use clap::{ArgGroup, Parser};
use rug::Integer;

#[derive(Parser, Debug)]
#[command(version, about = "Lenstra (2002) variant of the AKS primality test")]
#[command(group(ArgGroup::new("range").args(["theorem", "phi", "lambda"])))]
struct Cli {
    /// The number to test (decimal). Read from standard input if omitted.
    n: Option<String>,

    /// 1 <= a <= floor(sqrt(phi(r)) log2 n): the theorem's bound (default).
    #[arg(long)]
    theorem: bool,
    /// 1 <= a <= phi(r) - 1 (dissertation Algorithm 2.4.2; not certified).
    #[arg(long)]
    phi: bool,
    /// 1 <= a <= lambda(r) - 1 (Egan 2005; not certified).
    #[arg(long)]
    lambda: bool,

    /// Polynomial multiplier for the squarings.
    #[arg(long = "mul", default_value = "kronecker")]
    multiplier: Multiplier,

    /// Threads for the congruence tests (default: one per core, or
    /// RAYON_NUM_THREADS). 1 runs serially.
    #[arg(short = 'j', long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..))]
    threads: Option<u16>,

    /// Append a CSV row (the C++ log format) to this file.
    #[arg(long, value_name = "FILE")]
    log: Option<PathBuf>,

    /// Print only the verdict line.
    #[arg(short, long)]
    quiet: bool,
}

fn read_n(arg: Option<String>) -> Result<Integer, String> {
    let text = match arg {
        Some(s) => s,
        None => {
            let mut s = String::new();
            io::stdin()
                .read_to_string(&mut s)
                .map_err(|e| format!("reading n: {e}"))?;
            s
        }
    };
    let text = text.trim();
    text.parse::<Integer>()
        .map_err(|_| format!("not an integer: {text:?}"))
}

/// The C++ log's Comments column, so old and new logs can be compared.
fn note(t: &Trace) -> String {
    match &t.reason {
        Reason::Small => format!("{} is prime", t.n),
        Reason::PerfectPower { .. } => format!("{} is a perfect power", t.n),
        Reason::Divisor(d) => format!("{d} is a divisor"),
        Reason::NoRBelowN => "no r < n with ord_r(n) > log2(n)^2; all a < n coprime to n".into(),
        Reason::CongruencesHold | Reason::CongruenceFails(_) => {
            let p = t.params.expect("params");
            let end = match t.reason {
                Reason::CongruenceFails(a) => a,
                _ => 0,
            };
            format!(
                "a = {} ({}); End: a = {end}; r = {}; phi(r) = {}; lambda(r) = {}",
                t.a_max.expect("a_max"),
                t.options.range,
                p.r,
                p.phi,
                p.lambda
            )
        }
    }
}

fn append_log(path: &Path, t: &Trace) -> io::Result<()> {
    let new = !path.exists() || path.metadata()?.len() == 0;
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let mut f = OpenOptions::new().create(true).append(true).open(path)?;
    if new {
        writeln!(f, "Int, Cores, Prime (T/F), Time (milliseconds), Comments")?;
    }
    writeln!(
        f,
        "{}, {}, {}, {}, {}",
        t.n,
        t.threads,
        if t.prime { "T" } else { "F" },
        t.elapsed.as_millis(),
        note(t)
    )
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let n = match read_n(cli.n) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("aks: {e}");
            return ExitCode::from(2);
        }
    };
    let range = if cli.phi {
        ARange::PhiMinusOne
    } else if cli.lambda {
        ARange::LambdaMinusOne
    } else {
        ARange::Theorem
    };
    let options = Options {
        range,
        multiplier: cli.multiplier,
        threads: cli.threads.map(usize::from),
    };

    let t = match aks(&n, options) {
        Ok(t) => t,
        Err(e) => {
            println!("{e}.");
            return ExitCode::from(1);
        }
    };

    if !cli.quiet {
        if let Some(p) = t.params {
            println!("r = {}", p.r);
            println!("ord_r(n) = {}", p.order);
            println!("Phi({0}) = {1}, Lambda({0}) = {2}", p.r, p.phi, p.lambda);
            println!("Testing 1 <= a <= {} ({})", t.a_max.unwrap(), range);
            println!(
                "Multiplier: {}, threads: {}",
                options.multiplier.name(),
                t.threads
            );
        }
        match &t.reason {
            Reason::PerfectPower { base, exp } => {
                println!("{} = {base}^{exp} is a perfect power.", t.n)
            }
            Reason::Divisor(d) => println!("{d} is a divisor."),
            Reason::CongruenceFails(a) => println!("The a which fails is {a}"),
            _ => {}
        }
    }
    let uncertified = if t.prime && !t.certified() {
        " (not certified)"
    } else {
        ""
    };
    println!(
        "{} is {}{uncertified}.",
        t.n,
        if t.prime { "prime" } else { "not prime" }
    );
    if !cli.quiet {
        println!("Time taken: {} milliseconds", t.elapsed.as_millis());
    }

    if let Some(path) = &cli.log {
        if let Err(e) = append_log(path, &t) {
            eprintln!("aks: writing log {}: {e}", path.display());
        }
    }
    if t.prime {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(1)
    }
}
