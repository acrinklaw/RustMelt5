extern crate rust_melt_5;

use clap::Parser;
use crate::rust_melt_5::compute_thermodynamics;

#[derive(Parser, Debug)]
#[command(author="Austin Crinklaw", version = "0.1.0", about = "Compute melting temperature of mRNA/RNA duplexes")]
struct Args {
    // Input sequence
    #[arg(short, long)]
    sequence: String,
    // Sodium concentration
    #[arg(short, long)]
    na_conc: f64,
    // dNTP concentration
    #[arg(short, long)]
    rna_conc: f64,
}

fn main() {
    let args = Args::parse();
    let tm = compute_thermodynamics(&args.sequence, args.rna_conc, args.na_conc);
    println!("Melting temperature: {:.2}°C", tm);
}
