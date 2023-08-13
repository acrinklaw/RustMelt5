use clap::Parser;

pub fn compute_melting_temperature(nucleotide_concentration: f64, enthalpy: f64, entropy: f64) -> f64 {
    let tm = enthalpy / (entropy + 1.99 * (nucleotide_concentration/ 4.0).ln()) - 273.15;
    tm
}

pub fn count_terminal_au(sequence: &str) -> u8 {
    if sequence.is_empty() {
        return 0;
    }

    let first_char = sequence.chars().next().unwrap_or_default();
    let last_char = sequence.chars().last().unwrap_or_default();

    let is_first_au = matches!(first_char, 'A' | 'U');
    let is_last_au = matches!(last_char, 'A' | 'U');

    let mut count = 0;

    if is_first_au { count += 1; }
    if is_last_au { count += 1; }

    count
}

pub fn calculate_terminal_au_penalty(sequence: &str) -> (f64, f64) {
    let count = count_terminal_au(sequence);
    let enthalpy_penalty: f64 = 3140.0;
    let entropy_penalty: f64 = 9.1;
    (enthalpy_penalty * count as f64, entropy_penalty * count as f64)
}

pub fn calculate_delta_s(na: f64, duplex_length: usize) -> f64 {
    let square: f64 = na.ln() * na.ln();
    let a: f64 = -0.075 * na.ln() + 0.012 * square;
    let b: f64 = 0.018 * square;
    let g: f64 = a + b / duplex_length as f64;
    g
}

pub fn get_enthalpy(pair: &str) -> f64 {
    match pair {
        "AA" => -7480.0,
        "AC" => -6320.0,
        "AG" => -13940.0,
        "AU" => -6330.0,
        "UU" => -5430.0,
        "UG" => -12140.0,
        "UC" => -9650.0,
        "UA" => -6470.0,
        "CC" => -8880.0,
        "CA" => -5210.0,
        "CG" => -9470.0,
        "CU" => -9590.0,
        "GA" => -5770.0,
        "GG" => -9660.0,
        "GC" => -11900.0,
        "GU" => -6620.0,
        _ => panic!("Invalid pair {}", pair),
    }

}

pub fn get_entropy(pair: &str) -> f64{
    match pair {
        "AA" => -22.3,
        "AC" => -15.2,
        "AG" => -39.1,
        "AU" => -17.7,
        "UU" => -14.5,
        "UG" => -32.9,
        "UC" => -25.0,
        "UA" => -17.0,
        "CC" => -19.7,
        "CA" => -10.7,
        "CG" => -23.0,
        "CU" => -23.89,
        "GA" => -11.9,
        "GG" => -22.1,
        "GC" => -26.3,
        "GU" => -14.3,
        _ => panic!("Invalid pair {}", pair),
    }
}

pub fn compute_thermodynamics(sequence: &str, nucleotide_concentration: f64, na_concentration: f64) -> f64{
    //initial values for turner2006 mrna/rna
    let mut enthalpy: f64 = 0.0;
    let mut entropy: f64 = 0.0;
    for i in 0..sequence.len() - 1 {
        let pair = &sequence[i..i+2];
        enthalpy += get_enthalpy(pair);
        entropy += get_entropy(pair);
    }
    entropy = compute_entropy_1mna(entropy, 0, sequence.len() as i32 - 1);
    let entropy_correction: f64 = -3.22 * ((sequence.len() as f64 - 1.0) * calculate_delta_s(na_concentration, sequence.len()));
    entropy += entropy_correction;
    let (enthalpy_penalty, entropy_penalty) = calculate_terminal_au_penalty(sequence);
    enthalpy += enthalpy_penalty;
    entropy += entropy_penalty;
    // -12800.0, -52.0 are the values for duplex initiation
    let tm = compute_melting_temperature(nucleotide_concentration, enthalpy - 12800.0 , entropy - 52.0);
    tm
}

pub fn compute_entropy_1mna(entropy: f64, pos1: i32, pos2: i32) -> f64 {
    // Returns the entropy value in 1M [Na+] from initial entropy value in 0.1M Na using sodium correction
    // from Santa Lucia 2004
    let correction_constant: f64 = 0.1;
    let entropy_1mna: f64 = entropy - 0.368 * (pos2 - pos1).abs() as f64 * correction_constant.ln();
    entropy_1mna
}

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
