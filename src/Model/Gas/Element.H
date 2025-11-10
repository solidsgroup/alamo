#ifndef MODEL_GAS_ELEMENT_H
#define MODEL_GAS_ELEMENT_H

#include <cstdint>

// Get atomic mass for any element by index or symbol name
//    AtomicWeights("N")  // returns 14.0071, Nitrogen
//    AtomicWeights(1)    // return 4.0026022, Helium 

namespace Model {
namespace Gas {

enum class Element : int {
    H , He, Li, Be, B , C , N , O , F , Ne, Na, Mg, Al, Si, P , S , 
    Cl, Ar, K , Ca, Sc, Ti, V , Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, 
    As, Se, Br, Kr, Rb, Sr, Y , Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, 
    In, Sn, Sb, Te, I , Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, 
    Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W , Re, Os, Ir, Pt, Au, Hg, 
    Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U , Np, Pu, Am, Cm, 
    Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, 
    Nh, Fl, Mc, Lv, Ts, Og, COUNT };

constexpr double _AtomicWeights[static_cast<int>(Element::COUNT)] = {
    1.00802				,//	H		Hydrogen
    4.0026022			,//	He	Helium
    6.946					,//	Li	Lithium
    9.01218315		,//	Be	Beryllium
    10.812				,//	B		Boron
    12.0112				,//	C		Carbon
    14.0071				,//	N		Nitrogen
    15.9991				,//	O		Oxygen
    18.9984031625	,//	F		Fluorine
    20.17976			,//	Ne	Neon
    22.989769282	,//	Na	Sodium
    24.3052				,//	Mg	Magnesium
    26.98153843		,//	Al	Aluminium
    28.0851				,//	Si	Silicon
    30.9737619985	,//	P		Phosphorus
    32.062				,//	S		Sulfur
    35.451				,//	Cl	Chlorine
    39.9516				,//	Ar	Argon
    39.09831			,//	K		Potassium
    40.0784				,//	Ca	Calcium
    44.9559074		,//	Sc	Scandium
    47.8671				,//	Ti	Titanium
    50.94151			,//	V		Vanadium
    51.99616			,//	Cr	Chromium
    54.9380432		,//	Mn	Manganese
    55.8452				,//	Fe	Iron
    58.9331943		,//	Co	Cobalt
    58.69344			,//	Ni	Nickel
    63.5463				,//	Cu	Copper
    65.382				,//	Zn	Zinc
    69.7231				,//	Ga	Gallium
    72.6308				,//	Ge	Germanium
    74.9215956		,//	As	Arsenic
    78.9718				,//	Se	Selenium
    79.9043				,//	Br	Bromine
    83.7982				,//	Kr	Krypton
    85.46783			,//	Rb	Rubidium
    87.621				,//	Sr	Strontium
    88.9058382		,//	Y		Yttrium
    91.2223				,//	Zr	Zirconium
    92.906371			,//	Nb	Niobium
    95.951				,//	Mo	Molybdenum
    97.0					,//  Tc	Technetium
    101.072				,//	Ru	Ruthenium
    102.905492		,//	Rh	Rhodium
    106.421				,//	Pd	Palladium
    107.86822			,//	Ag	Silver
    112.4144			,//	Cd	Cadmium
    114.8181			,//	In	Indium
    118.7107			,//	Sn	Tin
    121.7601			,//	Sb	Antimony
    127.603				,//	Te	Tellurium
    126.904473		,//	I		Iodine
    131.2936			,//	Xe	Xenon
    132.905451966	,//	Cs	Caesium
    137.3277			,//	Ba	Barium
    138.905477		,//	La	Lanthanum
    140.1161			,//	Ce	Cerium
    140.907661		,//	Pr	Praseodymium  
    144.2423			,//	Nd	Neodymium
    145.0					,//  Pm	Promethium
    150.362				,//	Sm	Samarium
    151.9641			,//	Eu	Europium
    157.2492			,//	Gd	Gadolinium
    158.9253547		,//	Tb	Terbium
    162.5001			,//	Dy	Dysprosium
    164.9303295		,//	Ho	Holmium
    167.2593			,//	Er	Erbium
    168.9342195		,//	Tm	Thulium
    173.0451			,//	Yb	Ytterbium
    174.966695		,//	Lu	Lutetium
    178.4866			,//	Hf	Hafnium
    180.947882		,//	Ta	Tantalum
    183.841				,//	W		Tungsten
    186.2071			,//	Re	Rhenium
    190.233				,//	Os	Osmium
    192.2172			,//	Ir	Iridium
    195.0849			,//	Pt	Platinum
    196.9665704		,//	Au	Gold
    200.5923			,//	Hg	Mercury
    204.381				,//	Tl	Thallium
    207.211		  	,//	Pb	Lead
    208.980401		,//	Bi	Bismuth
    209.0					,//  Po	Polonium
    210.0					,//  At	Astatine
    222.0					,//  Rn	Radon
    223.0					,//  Fr	Francium
    226.0					,//  Ra	Radium
    227.0					,//  Ac	Actinium
    232.03774			,//	Th	Thorium
    231.035881		,//	Pa	Protactinium
    238.028913		,//	U		Uranium
    237.0					,//  Np	Neptunium
    244.0					,//  Pu	Plutonium
    243.0					,//  Am	Americium
    247.0					,//  Cm	Curium
    247.0					,//  Bk	Berkelium
    251.0					,//  Cf	Californium
    252.0					,//  Es	Einsteinium
    257.0					,//  Fm	Fermium
    258.0					,//  Md	Mendelevium
    259.0					,//  No	Nobelium
    262.0					,//  Lr	Lawrencium
    267.0					,//  Rf	Rutherfordium
    270.0					,//  Db	Dubnium
    269.0					,//  Sg	Seaborgium
    270.0					,//  Bh	Bohrium
    270.0					,//  Hs	Hassium
    278.0					,//  Mt	Meitnerium
    281.0					,//  Ds	Darmstadtium
    281.0					,//  Rg	Roentgenium
    285.0					,//  Cn	Copernicium
    286.0					,//  Nh	Nihonium
    289.0					,//  Fl	Flerovium
    289.0					,//  Mc	Moscovium
    293.0					,//  Lv	Livermorium
    293.0					,//  Ts	Tennessine
    294.0					 //  Og	Oganesson
};

constexpr unsigned int hash(const char* str) {
    unsigned int h = 0;
    while (*str) h = h * 31 + static_cast<unsigned char>(*str++);
    return h;
}

constexpr Element symbolToElement(const char* s) {
    switch (hash(s)) {
        case hash("H"):   return Element::H;
        case hash("He"):  return Element::He;
        case hash("Li"):  return Element::Li;
        case hash("Be"):  return Element::Be;
        case hash("B"):   return Element::B;
        case hash("C"):   return Element::C;
        case hash("N"):   return Element::N;
        case hash("O"):   return Element::O;
        case hash("F"):   return Element::F;
        case hash("Ne"):  return Element::Ne;
        case hash("Na"):  return Element::Na;
        case hash("Mg"):  return Element::Mg;
        case hash("Al"):  return Element::Al;
        case hash("Si"):  return Element::Si;
        case hash("P"):   return Element::P;
        case hash("S"):   return Element::S;
        case hash("Cl"):  return Element::Cl;
        case hash("Ar"):  return Element::Ar;
        case hash("K"):   return Element::K;
        case hash("Ca"):  return Element::Ca;
        case hash("Sc"):  return Element::Sc;
        case hash("Ti"):  return Element::Ti;
        case hash("V"):   return Element::V;
        case hash("Cr"):  return Element::Cr;
        case hash("Mn"):  return Element::Mn;
        case hash("Fe"):  return Element::Fe;
        case hash("Co"):  return Element::Co;
        case hash("Ni"):  return Element::Ni;
        case hash("Cu"):  return Element::Cu;
        case hash("Zn"):  return Element::Zn;
        case hash("Ga"):  return Element::Ga;
        case hash("Ge"):  return Element::Ge;
        case hash("As"):  return Element::As;
        case hash("Se"):  return Element::Se;
        case hash("Br"):  return Element::Br;
        case hash("Kr"):  return Element::Kr;
        case hash("Rb"):  return Element::Rb;
        case hash("Sr"):  return Element::Sr;
        case hash("Y"):   return Element::Y;
        case hash("Zr"):  return Element::Zr;
        case hash("Nb"):  return Element::Nb;
        case hash("Mo"):  return Element::Mo;
        case hash("Tc"):  return Element::Tc;
        case hash("Ru"):  return Element::Ru;
        case hash("Rh"):  return Element::Rh;
        case hash("Pd"):  return Element::Pd;
        case hash("Ag"):  return Element::Ag;
        case hash("Cd"):  return Element::Cd;
        case hash("In"):  return Element::In;
        case hash("Sn"):  return Element::Sn;
        case hash("Sb"):  return Element::Sb;
        case hash("Te"):  return Element::Te;
        case hash("I"):   return Element::I;
        case hash("Xe"):  return Element::Xe;
        case hash("Cs"):  return Element::Cs;
        case hash("Ba"):  return Element::Ba;
        case hash("La"):  return Element::La;
        case hash("Ce"):  return Element::Ce;
        case hash("Pr"):  return Element::Pr;
        case hash("Nd"):  return Element::Nd;
        case hash("Pm"):  return Element::Pm;
        case hash("Sm"):  return Element::Sm;
        case hash("Eu"):  return Element::Eu;
        case hash("Gd"):  return Element::Gd;
        case hash("Tb"):  return Element::Tb;
        case hash("Dy"):  return Element::Dy;
        case hash("Ho"):  return Element::Ho;
        case hash("Er"):  return Element::Er;
        case hash("Tm"):  return Element::Tm;
        case hash("Yb"):  return Element::Yb;
        case hash("Lu"):  return Element::Lu;
        case hash("Hf"):  return Element::Hf;
        case hash("Ta"):  return Element::Ta;
        case hash("W"):   return Element::W;
        case hash("Re"):  return Element::Re;
        case hash("Os"):  return Element::Os;
        case hash("Ir"):  return Element::Ir;
        case hash("Pt"):  return Element::Pt;
        case hash("Au"):  return Element::Au;
        case hash("Hg"):  return Element::Hg;
        case hash("Tl"):  return Element::Tl;
        case hash("Pb"):  return Element::Pb;
        case hash("Bi"):  return Element::Bi;
        case hash("Po"):  return Element::Po;
        case hash("At"):  return Element::At;
        case hash("Rn"):  return Element::Rn;
        case hash("Fr"):  return Element::Fr;
        case hash("Ra"):  return Element::Ra;
        case hash("Ac"):  return Element::Ac;
        case hash("Th"):  return Element::Th;
        case hash("Pa"):  return Element::Pa;
        case hash("U"):   return Element::U;
        case hash("Np"):  return Element::Np;
        case hash("Pu"):  return Element::Pu;
        case hash("Am"):  return Element::Am;
        case hash("Cm"):  return Element::Cm;
        case hash("Bk"):  return Element::Bk;
        case hash("Cf"):  return Element::Cf;
        case hash("Es"):  return Element::Es;
        case hash("Fm"):  return Element::Fm;
        case hash("Md"):  return Element::Md;
        case hash("No"):  return Element::No;
        case hash("Lr"):  return Element::Lr;
        case hash("Rf"):  return Element::Rf;
        case hash("Db"):  return Element::Db;
        case hash("Sg"):  return Element::Sg;
        case hash("Bh"):  return Element::Bh;
        case hash("Hs"):  return Element::Hs;
        case hash("Mt"):  return Element::Mt;
        case hash("Ds"):  return Element::Ds;
        case hash("Rg"):  return Element::Rg;
        case hash("Cn"):  return Element::Cn;
        case hash("Nh"):  return Element::Nh;
        case hash("Fl"):  return Element::Fl;
        case hash("Mc"):  return Element::Mc;
        case hash("Lv"):  return Element::Lv;
        case hash("Ts"):  return Element::Ts;
        case hash("Og"):  return Element::Og;
        default: return Element::H;
    }
}

inline double AtomicWeights(const char* symbol) {
    Element e = symbolToElement(symbol);
    return _AtomicWeights[static_cast<int>(e)];
}

inline double AtomicWeights(const int symbol) {
    return _AtomicWeights[symbol];
}

} // namespace Gas
} // namespace Model
#endif
