pub mod common;
mod methods;
use crate::common::{NucleotidePair, Method};
use crate::methods::NearestNeighbor;


pub fn nucleotide_pair_from_str(pair: &str) -> NucleotidePair {
    match pair {
        "AA" => NucleotidePair::AA,
        "AC" => NucleotidePair::AC,
        "AG" => NucleotidePair::AG,
        "AU" => NucleotidePair::AU,
        "UU" => NucleotidePair::UU,
        "UG" => NucleotidePair::UG,
        "UC" => NucleotidePair::UC,
        "UA" => NucleotidePair::UA,
        "CC" => NucleotidePair::CC,
        "CA" => NucleotidePair::CA,
        "CG" => NucleotidePair::CG,
        "CU" => NucleotidePair::CU,
        "GA" => NucleotidePair::GA,
        "GG" => NucleotidePair::GG,
        "GC" => NucleotidePair::GC,
        "GU" => NucleotidePair::GU,
        _ => panic!("Invalid pair {}", pair),
    }
}

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

pub fn compute_thermodynamics(sequence: &str, nucleotide_concentration: f64, na_concentration: f64) -> f64{
    //initial values for turner2006 mrna/rna
    let mut enthalpy: f64 = 0.0;
    let mut entropy: f64 = 0.0;

    // TODO: replace NearestNeighbor with method of choice
    let method = NearestNeighbor;
    for i in 0..sequence.len() - 1 {
        let pair_str = &sequence[i..i+2];
        let pair = nucleotide_pair_from_str(pair_str);

        let thermo_props = method.from_pair(&pair);
        enthalpy += thermo_props.enthalpy;
        entropy += thermo_props.entropy;
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
