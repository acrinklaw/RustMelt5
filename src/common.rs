pub enum NucleotidePair {
    AA, AC, AG, AU, UU, UG, UC, UA, CC, CA, CG, CU, GA, GG, GC, GU
}

pub struct ThermoProperties {
    pub enthalpy: f64,
    pub entropy: f64,
}

pub trait Method {
    fn from_pair(&self, pair: &NucleotidePair) -> ThermoProperties;
}