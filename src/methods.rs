use crate::common::{NucleotidePair, ThermoProperties, Method};

pub struct NearestNeighbor;

impl Method for NearestNeighbor {
    fn from_pair(&self, pair: &NucleotidePair) -> ThermoProperties {
        match pair {
            NucleotidePair::AA => ThermoProperties { enthalpy: -7480.0, entropy: -22.3 },
            NucleotidePair::AC => ThermoProperties { enthalpy: -6320.0, entropy: -15.2 },
            NucleotidePair::AG => ThermoProperties { enthalpy: -13940.0, entropy: -39.1 },
            NucleotidePair::AU => ThermoProperties { enthalpy: -6330.0, entropy: -17.7 },
            NucleotidePair::UU => ThermoProperties { enthalpy: -5430.0, entropy: -14.5 },
            NucleotidePair::UG => ThermoProperties { enthalpy: -12140.0, entropy: -32.9 },
            NucleotidePair::UC => ThermoProperties { enthalpy: -9650.0, entropy: -25.0 },
            NucleotidePair::UA => ThermoProperties { enthalpy: -6470.0, entropy: -17.0 },
            NucleotidePair::CC => ThermoProperties { enthalpy: -8880.0, entropy: -19.7 },
            NucleotidePair::CA => ThermoProperties { enthalpy: -5210.0, entropy: -10.7 },
            NucleotidePair::CG => ThermoProperties { enthalpy: -9470.0, entropy: -23.0 },
            NucleotidePair::CU => ThermoProperties { enthalpy: -9590.0, entropy: -23.89 },
            NucleotidePair::GA => ThermoProperties { enthalpy: -5770.0, entropy: -11.9 },
            NucleotidePair::GG => ThermoProperties { enthalpy: -9660.0, entropy: -22.1 },
            NucleotidePair::GC => ThermoProperties { enthalpy: -11900.0, entropy: -26.3 },
            NucleotidePair::GU => ThermoProperties { enthalpy: -6620.0, entropy: -14.3 },
        }
    }
}