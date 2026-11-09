//! CLI binary `rho-factor`.
//!
//! Usage: `rho-factor <N> [--threads T] [--batch-size B]`
//!
//! Factors an integer N below 2^80 using parallel Pollard rho (Brent + batched GCD).
//! Prints the smaller factor and its cofactor on success, or `N is prime` (exit 0) when N is
//! certified prime.  Invalid or unsupported input and an exhausted budget are reported
//! separately on stderr (exit 1): exhaustion says nothing about whether N is prime.

use clap::Parser;
use rho::factor::{FactorError, PRIME_PROOF_BOUND};

#[derive(Parser, Debug)]
#[command(name = "rho-factor", about = "Pollard rho integer factorization")]
struct Args {
    /// Integer to factor (decimal, below 2^80).
    n: u128,

    /// Number of c-values to try in parallel (defaults to logical CPU count).
    #[arg(long, default_value_t = 0)]
    threads: u64,

    /// GCD batch size for Montgomery's batched-GCD trick (default: 128).
    #[arg(long, default_value_t = 128)]
    batch_size: usize,
}

fn main() {
    let args = Args::parse();
    let n = args.n;

    // Default threads: number of available logical CPUs.
    let max_c = if args.threads == 0 {
        rayon::current_num_threads() as u64
    } else {
        args.threads
    };

    match rho::factor::factor(n, max_c, args.batch_size) {
        Ok(found) => {
            let p = found.d;
            let q = n / p;
            let (lo, hi) = if p <= q { (p, q) } else { (q, p) };
            println!("{n} = {lo} × {hi}");
        }
        Err(FactorError::Prime) => println!("{n} is prime"),
        Err(FactorError::InvalidInput) => {
            eprintln!("error: N must be > 1");
            std::process::exit(1);
        }
        Err(FactorError::Unsupported) => {
            eprintln!(
                "rho-factor: {n} is outside the supported range (N < 2^80, and an odd N above \
                 {PRIME_PROOF_BOUND} that tests prime cannot be certified)"
            );
            std::process::exit(1);
        }
        Err(FactorError::BudgetExhausted { constants }) => {
            eprintln!(
                "rho-factor: no factor of {n} found after {constants} c-values; this does not \
                 show {n} is prime (try --threads with a larger value)"
            );
            std::process::exit(1);
        }
    }
}
