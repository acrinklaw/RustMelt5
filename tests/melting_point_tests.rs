extern crate rust_melt_5;
use crate::rust_melt_5::compute_thermodynamics;

#[test]
fn test_compute_thermodynamics_no_terminal_au() {
    let test_sequence = "GAGUCAAACAAGGAC"; 
    let test_nucleotide_concentration = 2e-9; 
    let test_na_concentration = 0.15;

    let result = compute_thermodynamics(test_sequence, test_nucleotide_concentration, test_na_concentration);
    let formatted_result = format!("{:.2}", result).parse::<f64>().unwrap();
    let expected_output = 57.99; 
    assert_eq!(formatted_result, expected_output);
}

#[test]
fn test_compute_thermodynamics_with_one_terminal_au() {
    let test_sequence = "GAGUCAAACAAGGAA"; 
    let test_nucleotide_concentration = 2e-9; 
    let test_na_concentration = 0.15;

    let result = compute_thermodynamics(test_sequence, test_nucleotide_concentration, test_na_concentration);
    let formatted_result = format!("{:.2}", result).parse::<f64>().unwrap();
    let expected_output = 54.44; 
    assert_eq!(formatted_result, expected_output);
}

#[test]
fn test_compute_thermodynamics_with_two_terminal_au() {
    let test_sequence = "AAGUCAAACAAGGAA"; 
    let test_nucleotide_concentration = 2e-9; 
    let test_na_concentration = 0.15;

    let result = compute_thermodynamics(test_sequence, test_nucleotide_concentration, test_na_concentration);
    let formatted_result = format!("{:.2}", result).parse::<f64>().unwrap();
    let expected_output = 49.46; 
    assert_eq!(formatted_result, expected_output);
}

