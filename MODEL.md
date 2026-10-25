# The Reef Twin simulation model

*Last updated 2026-10-24 (species catalog 2026.10.7).*

Reef Twin is a reef aquarium simulator. This document is the model behind it: the equations, every constant, where each number came from, and the tests the model has to pass. It is published so that anyone can check the work and tell me where it is wrong.

## How to read this

**What the model is for.** It is meant to get the direction and the rough size of things right: how fast alkalinity falls with a given coral load, how long a new tank takes to cycle, what house CO2 does to pH, how a heat wave or an overdose plays out. It is not a prediction of your own tank to two decimal places. No two real tanks are the same, and many of the numbers below are honest estimates.

**Not every part is equally certain.** Seawater chemistry is decades-old oceanography and is checked against the standard tools oceanographers use. Bacteria, coral biology and pests are much less certain. Every number carries one of four tags so you can see which is which:

| Tag | Meaning |
|---|---|
| **S** | From a peer-reviewed or primary source, listed in §2.14 |
| **H** | Hobby consensus with a named reference (Holmes-Farley, Sprung & Delbeek, Borneman, BRS tests, well-known forum threads) |
| **E** | An estimate chosen to reproduce known behaviour; no single citable number exists |
| **T** | A tuning value with no physical meaning, chosen so a validation scenario passes |

§2.13 lists the claims that could not be fully sourced. If you can source one of them, or think it is wrong, that is the most useful thing you can tell me.

**Validation.** §2.12 lists 18 scenarios the model runs automatically on every change, such as a new-tank cycle, house CO2, a kalkwasser snowstorm, low magnesium, a heat wave and a power outage. All 18 pass today.

**How to challenge it.** Open an issue on this repository, or reply in the forum thread where you found this. Say what the model does, what you would expect, and your source or experience. Changes made because of feedback are dated in the text.

**Reading tips.** Units are shown the way the hobby uses them (°F, gallons, dKH, ppm), though the maths is done in metric. "§" numbers refer to sections of this document. You can skip the equations and read the tables and the "Check" and "Required behavior" paragraphs: they compare the model with what reefers actually see.

---

## Contents

- [2.1 Conventions](#21-conventions)
- [2.2 State vector](#22-state-vector): what the model keeps track of
- [2.3 Temperature, evaporation and salinity](#23-temperature-evaporation-and-salinity)
- [2.4 Carbonate chemistry](#24-carbonate-chemistry): alkalinity, pH, calcium, magnesium
- [2.5 Nitrogen and phosphorus](#25-nitrogen-and-phosphorus)
- [2.6 Light](#26-light)
- [2.7 Coral biology](#27-coral-biology)
- [2.8 Fish biology](#28-fish-biology)
- [2.8a Coralline algae](#28a-coralline-algae-rock-surface-as-built-2026-10-19)
- [2.8b Cleanup crew and display copper](#28b-cleanup-crew-and-display-copper-as-built-2026-10-23)
- [2.9 Nuisances and events](#29-nuisances-and-events-optional-extensions-off-by-default)
- [2.10 Equipment models](#210-equipment-models)
- [2.11 Species parameter schema](#211-species-parameter-schema)
- [2.12 Validation scenarios](#212-validation-scenarios)
- [2.13 Claims that could not be fully sourced](#213-claims-that-could-not-be-fully-sourced)
- [2.14 References](#214-references)

---

## 2.1 Conventions

- Internal units are SI: °C, kg, L, mol. The UI shows °F, gallons, dKH, ppm, as the hobby does.
- Concentration display follows hobby test kits: nitrate as NO3 ppm (not NO3-N; ×4.43), phosphate as PO4 ppm (not P; ×3.07), ammonia as total ammonia nitrogen (TAN) in mg/L N shown as "ammonia (NH3/NH4) ppm", alkalinity in dKH (1 dKH = 0.3571 meq/L = 17.86 ppm CaCO3), pH on the NBS scale because that is what hobby probes and kits read.
- Time is simulated days. One integration step is at most 0.02 d (28.8 min). A step is subdivided whenever any state variable would change by more than 2 % of its display range in one step, or the carbonate solver fails to converge.
- The model is a pure function `advance(state, catalog, schedule, dt, rngSeed)`. All randomness is drawn from a counter-based generator keyed on `(tankId, stepIndex)` so advancing 30 days in one call or 30 calls gives the same result. That is what lets a tank nobody is watching be brought up to date in one jump.

## 2.2 State vector

### Tank and water (one per aquarium)

| Symbol | Meaning | Unit | Initial (new saltwater) |
|---|---|---|---|
| V | display + sump water volume | L | from dimensions, minus rock and sand displacement |
| A_s | free water surface area (display + sump) | m² | from dimensions |
| A_g | glass area in contact with room air | m² | from dimensions |
| T | water temperature | °C | 25.5 |
| S | salinity | psu | 35 |
| TA | total alkalinity | µeq/kg | 2,860 (8.0 dKH) |
| DIC | dissolved inorganic carbon | µmol/kg | solved from TA and air pCO2 at start |
| Ca | calcium | mg/L | 420 |
| Mg | magnesium | mg/L | 1,300 |
| TAN | total ammonia nitrogen | mg N/L | 0 |
| NO2 | nitrite nitrogen | mg N/L | 0 |
| NO3 | nitrate nitrogen | mg N/L | 0 |
| PO4 | orthophosphate phosphorus | mg P/L | 0.005 |
| DOC | labile organic carbon (food, waste, mucus) | mg C/L | 0 |
| P_sorb | phosphate bound to rock and sand surfaces | mg P | 0 |
| O2 | dissolved oxygen | mg/L | saturation |
| B_AOB, B_AOA, B_NOB | ammonia-oxidising bacteria and archaea and nitrite-oxidiser biomass as fraction of surface carrying capacity | 0..1 | 10⁻⁵, 10⁻⁶, 10⁻⁶ on dry rock (B_AOB = B_NOB = 0.05 with bottled bacteria) |
| B_het | heterotroph biomass index for carbon dosing | 0..1 | 0.01 |
| D_cap | denitrification capacity index (matures with rock and sand age) | 0..1 | 0 |
| Alg, Cyano, Dino, Aip | nuisance biomass indices (optional extensions) | 0..1 | 0 |
| GFO_cap, GAC_cap | remaining media capacity | mg, mg | per media change |
| pH, Ω_arag, pCO2 | derived each step, not state | | |

### Per coral colony

| Symbol | Meaning | Unit |
|---|---|---|
| L | characteristic size (branching: max dimension; plating: diameter; LPS: head count × head size) | cm |
| M | skeletal CaCO3 mass | g |
| H | health | 0..1 |
| Z | zooxanthellae density, normalized to species baseline | ~1 |
| C | host pigment (GFP/chromoprotein) expression index | 0..1 |
| I_acc | acclimated irradiance (14-day mean daily PAR) | µmol m⁻² s⁻¹ |
| DD | thermal stress accumulator | °C·d |
| age, maturity, dead, skeletonAge | | |

### Per fish

| Symbol | Meaning | Unit |
|---|---|---|
| L | standard length | cm |
| W | mass, from length-weight relation W = a L^b (FishBase) | g |
| H | health | 0..1 |
| sat | satiation | 0..1 |
| age | | d |
| qt | days in quarantine before introduction | d |
| infected | pathogen load index (ich, velvet) | 0..1 |

## 2.3 Temperature, evaporation and salinity

### Energy balance

ρ V c_p dT/dt = P_heater + P_pumps + η_L P_lights − U A_g (T − T_room) − h_s A_s (T − T_air) − h_fg E − P_chiller

| Constant | Value | Tag | Source |
|---|---|---|---|
| ρ seawater | 1,025 kg/m³ | S | Sharqawy, Lienhard & Zubair 2010 |
| c_p seawater | 3,990 J/(kg·K) | S | Sharqawy et al. 2010 |
| U glass-to-room | 6 W/(m²·K) for 10–12 mm glass, still air | S (derived) | Series of water-side convection (~300), glass conduction (1.0 W/(m·K) / 0.012 m ≈ 83), outside natural convection (~7 W/(m²·K)); Incropera & DeWitt |
| h_s surface convection | 8 W/(m²·K) still, 15 with fan | E | Natural convection range, Incropera & DeWitt |
| h_fg latent heat | 2.43 MJ/kg | S | Standard |
| η_L light heat fraction into water | LED 0.35, T5 0.5, metal halide 0.6 | E | Fraction of electrical input ending as absorbed radiant energy in the water; LED driver and heatsink losses go to air |
| Pump heat fraction | submerged 0.95, external 0.5 | E | Nearly all shaft work becomes heat in the water |

Heater: thermostat with ±0.3 °C hysteresis around setpoint, output = rated watts when on (H: typical hobby heater behavior). Chiller: removes rated capacity (BTU/h × 0.293 = W) when T > setpoint + 0.5 °C until T < setpoint − 0.5 °C.

**Sanity check (120 gal, 2 × 200 W LED at 10 h, 100 W pumps, room 22 °C, tank 25.5 °C):** glass loss ≈ 6 × 2.2 × 3.5 ≈ 46 W; evaporative loss ≈ 6.9 kg/d × 2.43 MJ ≈ 194 W; lights average 58 W; pumps 95 W; heater must average ≈ 90 W, a 25 % duty cycle on a 300 W heater. This matches the hobby rule of 3 to 5 W per gallon of heater capacity (H: BRS heater guide).

### Evaporation

E [kg/h] = A_s · (25 + 19 v) · (x_s − x_a)

where v is air speed over the surface (m/s; 0.1 still room, 1.0 with a fan or open-top with strong flow), x_s is the saturation humidity ratio at water temperature reduced by 2 % for seawater, and x_a is the room humidity ratio from T_room and RH. S: ASHRAE/Carrier evaporation correlation as given in Engineering ToolBox "Evaporation from a water surface".

**Check:** 120 gal open top, 0.74 m² surface, water 26 °C, room 22 °C at 40 % RH, still air: E ≈ 0.39 kg/(m²·h) ≈ 6.9 L/day ≈ 1.5 % of volume per day. Hobby experience is 1 to 2 % per day for open tops (H: common Reef2Reef reports). A lid reduces E by 70 to 90 % (E).

### Salinity

Salt mass is conserved. Without top-off: dS/dt = S · E / (ρ V) and dV/dt = −E/ρ. With ATO: V held, S unchanged. Kalkwasser top-off: as ATO plus the dose in §2.4. Water change of fraction f with mix salinity S_m: S ← (1 − f) S + f S_m, and the same linear blend for every dissolved variable. Specific gravity for display: SG ≈ 1 + 0.000754 S (fits 35 psu → 1.0264 at 25/25 °C; H: Holmes-Farley "Specific gravity: oh how complicated"). Refractometer calibration error is an optional setting (±0.001) because so many real tanks run at the wrong salinity.

## 2.4 Carbonate chemistry

### Why TA and DIC are the state, not pH

pH is an outcome of alkalinity, dissolved CO2, temperature and salinity. Treating it as a state with ad-hoc adjustments (as early versions of this project did) cannot reproduce the four behaviors every experienced reefer knows: house air CO2 depresses pH, kalkwasser raises it, a calcium reactor lowers it, and photosynthesis swings it. All four fall out of the standard seawater CO2 system (S: Zeebe & Wolf-Gladrow 2001; Dickson, Sabine & Christian 2007 "Guide to best practices for ocean CO2 measurements"; Holmes-Farley "The relationship between alkalinity and pH" and "Low pH: causes and cures").

### Speciation

Given TA, DIC, T, S, solve for [H⁺] from

TA = [HCO3⁻] + 2[CO3²⁻] + [B(OH)4⁻] + [OH⁻] − [H⁺] (phosphate and silicate alkalinity added when PO4 > 0.1 mg/L)

with [HCO3⁻] = DIC · K1 h / (h² + K1 h + K1 K2), [CO3²⁻] = DIC · K1 K2 / (h² + K1 h + K1 K2), [B(OH)4⁻] = B_T K_B / (K_B + h), [OH⁻] = K_W / h.

| Constant | Source | Tag |
|---|---|---|
| K1, K2 | Lueker, Dickson & Keeling 2000 (total scale, S 19–43, T 2–35 °C) | S |
| K_B | Dickson 1990 | S |
| K_W | Millero 1995 | S |
| K_0 (CO2 solubility) | Weiss 1974 | S |
| B_T | 0.0004157 × S/35 mol/kg (Uppström 1974) | S |
| Ksp aragonite | Mucci 1983 | S |
| NBS vs total pH scale offset | +0.13 at S 35, 25 °C, applied for display | E (the offset is real and documented in Dickson et al. 2007; the single number is a simplification) |

Solve with Newton's method seeded by the previous step's h; two iterations suffice in practice. Then pCO2 = [CO2*]/K_0, Ω_arag = [Ca²⁺][CO3²⁻]/Ksp.

### Processes that change TA, DIC, Ca and Mg

| Process | Per mol | ΔTA (eq) | ΔDIC (mol) | ΔCa (mol) | ΔMg (mol) | Tag |
|---|---|---|---|---|---|---|
| Coral calcification | CaCO3 formed | −2 | −1 | −1 | −0.01 (aragonite is low-Mg) | S (stoichiometry); Mg fraction E |
| Coralline algae calcification | high-Mg calcite | −2 | −1 | −0.87 | −0.13 | S: coralline skeletons are 10–18 mol % Mg (Chave 1954); modeled as a fixed fraction of coral calcification, 15 % (E). Coverage is now modelled (§2.8a), but its calcification still uses this fixed term |
| Abiotic precipitation on heaters/pumps | CaCO3 | −2 | −1 | −1 | −0.05 | S (stoich.); rate below |
| Photosynthesis / respiration | C fixed or released | ±0.15 (nitrate assimilation) | −1 / +1 | 0 | 0 | S: Zeebe & Wolf-Gladrow |
| Nitrification NH4⁺ → NO3⁻ | N | −2 | 0 | 0 | 0 | S: standard |
| Ammonia excretion by fish (as NH3) | N | +1 | 0 | 0 | 0 | S |
| Denitrification NO3⁻ → N2 | N | +1 | +1 (from organic C) | 0 | 0 | S |
| Kalkwasser Ca(OH)2 | Ca | +2 | 0 | +1 | 0 | S; saturated limewater ≈ 0.0203 mol/L → 40.5 meq/L alk, 811 mg/L Ca (H: Holmes-Farley "What your grandmother never told you about lime") |
| Two-part, sodium bicarbonate | NaHCO3 | +1 | +1 | 0 | 0 | S |
| Two-part, sodium carbonate (soda ash) | Na2CO3 | +2 | +1 | 0 | 0 | S |
| Two-part, calcium chloride | CaCl2 | 0 | 0 | +1 | 0 | S |
| Magnesium (MgCl2 / MgSO4 mix) | Mg | 0 | 0 | 0 | +1 | S |
| Calcium reactor effluent | CaCO3 dissolved by CO2 | +2 | +2 | +1 | +0.05 (media is often high-Mg) | S: the extra DIC is why reactors depress pH (H: Holmes-Farley "Calcium reactors") |
| Carbon dosing (vinegar) | acetic acid respired | 0 net | +2 per acetate | 0 | 0 | S |
| Water change | linear blend to mix values | | | | | |

Dosing arithmetic shown to the user: ΔdKH = dose_mL × C_eq (meq/mL) / V (L) × 2.8; ΔCa (ppm) = dose_mL × C (mg/mL) / V (L). Catalog solutions carry their C_eq and C with a reference (BRS, Tropic Marin, Red Sea, ESV published strengths; H).

### Air-water CO2 exchange

F_CO2 [mol/d] = k_tot · A_s · K_0 · (pCO2_air − pCO2_water)

k_tot = k_surface + k_skimmer. k_surface = 0.02 m/h still surface to 0.10 m/h with strong surface agitation (S: order of magnitude from lake and tank gas-transfer studies, Wanninkhof 1992 at low wind; E for the mapping to flow setting). k_skimmer = ε Q_air / A_s with ε = 0.4 (E), Q_air the skimmer's air draw (catalog, L/h). Room pCO2 is a user setting with presets: outdoor 420 ppm, ventilated house 550, closed house in winter 900, basement fish room 1,200 (H: Holmes-Farley "Low pH" article discusses 1,000+ ppm indoors; E for the preset values). Running the skimmer air line outdoors sets pCO2_air for k_skimmer to outdoor.

**Required behavior (validation):** at 8 dKH, 25 °C, S 35: equilibrium pH 8.24 ± 0.05 NBS at 420 ppm, 7.93 ± 0.07 at 1,000 ppm (PyCO2SYS with the constants above gives 8.243 and 7.925; the NBS–total offset is 0.137); a skimmer-equipped 120 gal reaches 90 % of equilibrium within 6 to 12 h (H: Holmes-Farley's alkalinity/pH table; the equilibration time is E).

### Abiotic precipitation

R_p [mol CaCO3/d] = k_p · V · Σ max(0, Ω_eff − Ω_c)² summed over two mineral paths, with Ω_eff computed at heater-surface temperature (T + 8 °C on the element; E) for the heater's share (30 % of surfaces while the element is on; E) and at bulk T elsewhere. Aragonite: Ω_c,ar = 5.5 · (Mg/1300)^0.5 (E; Mg inhibits aragonite growth, S: Berner 1975; H: Holmes-Farley "Magnesium in reef aquaria"). Calcite: Ω_cal = 1.52 · Ω_ar (S: Mucci 1983 solubility ratio at 25 °C, S 35) against Ω_c,cal = 12 · (Mg/1300)^4 (E for the base, T for the exponent): magnesium poisons calcite nucleation so strongly that at normal Mg only a gross overdose precipitates, while at Mg 950 ordinary 9 dKH / Ca 450 water loses about 0.1 dKH/day on heaters and pump impellers, the low-magnesium signature (H). k_p is T, set so that Ca 500 / 12 dKH / pH 8.45 / Mg 1300 precipitates ≈ 0.5 dKH/day, the "snowstorm" regime hobbyists report after kalk overdoses (H); implemented value 1.25 × 10⁻⁶ mol L⁻¹ d⁻¹ per unit excess². Precipitation removes Ca and Mg at −1 and −0.05 per mol (E). Kalkwasser also co-precipitates phosphate on contact (§2.5 table).

## 2.5 Nitrogen and phosphorus

### Inputs

Feed F (g dry/day) with protein fraction PC and phosphorus fraction P_f (catalog per food type; typical pellet PC 0.45, P_f 0.012; H: label analyses):

- Ammonia produced by fish: P_TAN = 0.092 · F_eaten · PC [g N/day] (S: Timmons & Ebeling 2010, eq. for TAN production in recirculating aquaculture).
- Uneaten food and feces go to DOC: F_waste = F · (1 − e_eat) + F_eaten · 0.25 (E), mineralized to TAN at k_min = 0.15 /d (E) unless skimmed first.
- Phosphorus: P_in = F · P_f · 0.7 released as PO4 over 2 days (E).
- Other P inputs: RODI TDS creep (user setting, mg P/L of top-off water), rock leaching from P_sorb.

### Nitrification (Monod on surface-bound biomass)

r_AOB [g N/d] = V_max · C_surf · B_AOB · TAN/(K_TAN + TAN) · θ^(T−20) · O2/(K_O + O2)
r_AOA [g N/d] = V_max · C_surf · B_AOA · TAN/(K_TAN,AOA + TAN) · θ^(T−20) · O2/(K_O + O2)
r_NOB [g N/d] = V_max,NOB · C_surf · B_NOB · NO2/(K_NO2 + NO2) · θ^(T−20) · O2/(K_O + O2)

| Constant | Value | Tag | Source |
|---|---|---|---|
| µ_max AOB at 20 °C | 0.75 /d | S | Metcalf & Eddy (Tchobanoglous et al. 2014) kinetic coefficients for Nitrosomonas |
| µ_max NOB at 20 °C | 0.6 /d | E | NOB generally slower; Nitrospira dominate aquaria (S: Hovanec & DeLong 1996; Hovanec et al. 1998) |
| K_TAN | 0.7 mg N/L | S | Metcalf & Eddy |
| µ_max AOA at 20 °C | 0.32 /d | S (rate) / E (temperature) | 0.027 /h (0.65 /d) for *Nitrosopumilus maritimus* at 30 °C (Martens-Habbena et al. 2009), brought to 20 °C with θ |
| K_TAN,AOA | 0.002 mg N/L | S | 0.133 µM total ammonium for *N. maritimus*, 30 °C, pH 7.5 (Martens-Habbena et al. 2009); Jung et al. 2022 find lower total-ammonium Km at higher pH, so at reef pH this is conservative |
| K_NO2 | 0.2 mg N/L | S | Metcalf & Eddy |
| θ temperature coefficient | 1.072 | S | Metcalf & Eddy |
| K_O half-saturation for O2 | 0.5 mg/L | S | Metcalf & Eddy |
| Areal capacity V_max | 0.5 g TAN/(m²·d) at saturation | S (range) | Timmons & Ebeling: 0.1–0.9 g/(m²·d) for biofilter media |
| Effective nitrifying surface of live rock | 0.5 m²/kg | E | No reliable published figure for porous aquarium rock; chosen so 1 lb/gal rock is heavily over-provisioned once mature, as hobby experience says |
| Effective surface of sand | 20 m²/L of bed, 2 cm active depth | E | |
| Biomass growth | dB/dt = µ_max θ^(T−20) · B · (substrate term) · (1 − B) · ρ − b B, decay b = 0.05 /d for all three guilds | S (form) / E (decay) | Salem et al. 2006 measured 0.2 /d aerobic decay in activated sludge; biofilm decay is lower and unmeasured for aquaria |
| Seed biomass on dry rock | B_AOB 10⁻⁵, B_AOA 10⁻⁶, B_NOB 10⁻⁶; bottled bacteria B_AOB = B_NOB = 0.05, no archaea | E | Bottled products are bacterial; archaea reach new rock with water and livestock |

**Two ammonia-oxidising guilds** (added 2026-10-18 after tester feedback). A single Monod population sets the standing ammonium of a mature tank at TAN* = K·b / (µ θ^(T−20) f_O2 − b), whatever the load: more feeding grows more bacteria, not more ammonia. With the Nitrosomonas constants above that is ≈ 0.036 mg N/L at 25.5 °C, an order of magnitude above reef water (open ocean below 0.02 and typically 0.006 mg N/L; HIMB intake 0.0025 ± 0.0021; Holmes-Farley 2007), and tanks crept to it over months. Marine aquarium biofilters carry ammonia-oxidising archaea alongside AOB, in proportions that vary from tank to tank (archaea dominated five of eight hobby saltwater filters, Sauder et al. 2011; one coral tank swung between AOB- and AOA-dominated, Bagchi et al. 2014). The archaea's affinity is about 350 times higher, so they hold a mature tank near TAN* ≈ 0.0003 mg N/L; the faster AOB still carry a new tank's cycle, where ammonia is far above both K. Both guilds draw on the same ammonium; each step's oxidation is capped by what is available, and **ρ** is the share of the guilds' potential rate the available ammonium allowed in that step (1 unless capped), so biomass grows on what it actually oxidised. Without ρ the archaea's tiny K makes the step stiff and the biomass runs away. Nitrification's oxygen demand (4.57 g O2 per g N) is charged on the oxidised amount for the same reason. In the hourly readout a mature tank shows TAN of order 0.001 to 0.004 mg N/L, roughly in proportion to the feeding, because the sample falls just after a step's ammonia input. An antibiotic treatment (chemiclean) sets back AOB and NOB but not the archaea, since erythromycin targets bacterial ribosomes (E). Tank files saved before 2026-10-18 get B_AOA = B_AOB on load (E).

Toxic un-ionized ammonia fraction for stress: f_NH3 = 1 / (1 + 10^(pKa − pH)), pKa = 0.0901821 + 2729.92 / T_K (S: Emerson et al. 1975; seawater pKa is ~0.1 lower, E).

**Required behavior:** an uncycled 120 gal with 120 lb dry rock and 2 ppm TAN dosed shows a TAN peak, a nitrite peak lagging 7 to 14 d, and TAN and NO2 both below 0.1 mg/L between day 18 and day 35 (H: standard hobby cycling timelines; S: Hovanec & DeLong 1996 for the 3 to 5 week establishment). With bottled nitrifiers (B₀ = 0.1), under 10 d.

### Denitrification

r_den [g N/d] = V_den · D_cap · NO3/(K_den + NO3) · DOC/(K_C + DOC) · θ^(T−20), with V_den per kg rock plus per liter of sand below 5 cm depth, K_den = 1 mg N/L (E), K_C = 1 mg C/L (E), D_cap maturing from 0 to 1 over 6 months with τ = 60 d (H: "rock matures over months", Sprung & Delbeek). V_den is T, set so a mature 120 gal with 120 lb rock removes ≈ 0.7 ppm NO3/day at 10 ppm NO3 with moderate feeding, which is why such tanks sit at low single-digit nitrate without other export (H: common experience; S for the mechanism: Risk & Muller 1983 on nitrogen in porous reef rock; Rasheed et al. 2002 on carbonate sands).

### Carbon dosing

Dose of organic carbon C_d [mol C/d] (vinegar 5 %: 0.83 mol C/L; vodka 40 %: 6.9 mol C/L; S: stoichiometry). Heterotroph growth takes N and P at Redfield proportions: ΔNO3 = −C_d · (16/106) · y, ΔPO4 = −C_d · (1/106) · y, y = 0.5 growth yield (S: Redfield 1958 ratio; yield E). B_het rises with dosing; without a skimmer the biomass stays in the water (bacterial bloom event, O2 demand +); with a skimmer it is exported. Overdose (NO3 < 0.5 and PO4 < 0.01) triggers the "ultra-low nutrient" stress in §2.7 and raises dinoflagellate risk in §2.9.

As implemented (2026-10-18): uptake is limited by the scarcer of N and P (so a tank usually runs out of phosphate first, as hobbyists report), ammonium is taken before nitrate, carbon the bacteria cannot use is respired (one O2 per C), and without a skimmer the new biomass stays as organic carbon that decays back to ammonium; B_het is not tracked as its own index. Until then the code added the dose to DOC, whose mineralisation made ammonia, so carbon dosing raised nitrate (found while building test tanks).

### Phosphate sinks

| Sink | Model | Tag |
|---|---|---|
| Coral and algae uptake | proportional to calcification and algal growth at C:P 106:1 by carbon fixed | S (Redfield) |
| Skimmer | removes a fraction s of DOC per pass-equivalent; s such that 20–35 % of organic input is exported (S: Feldman et al. 2010–2011, Advanced Aquarist TOC studies) | S |
| GFO | ΔPO4 = −k_g · PO4 · GFO_remaining; capacity 25 mg PO4 per g GFO | H (BRS and Holmes-Farley "Iron oxide hydroxide (GFO) phosphate binders"); capacity E |
| Rock/sand sorption | reversible: dP_sorb/dt = k_a · PO4 · (P_cap − P_sorb) − k_d · P_sorb; explains leaching after a high-nutrient period | S (mechanism: Millero 2001 phosphate adsorption on calcite); constants E |
| Kalkwasser precipitation | ΔPO4 = −0.02 mg/L per 1 dKH added as kalk | H (Holmes-Farley on limewater reducing phosphate); number E |
| Lanthanum chloride | stoichiometric, LaPO4 | S |
| Refugium macroalgae | Chaetomorpha uptake 2 mg N and 0.13 mg P per g dry per day under 100 µmol light; growth doubling ≈ 10 d; harvest removes mass | E (order of magnitude from macroalgae nutrient uptake literature, e.g. Neori et al. 2004 on Ulva biofilters) |
| Water change | linear blend | |

### Oxygen

Tracked simply: dO2/dt = k_O2 (O2_sat(T,S) − O2) + photosynthesis − respiration(fish, corals, bacteria, DOC). O2_sat from Garcia & Gordon 1992 (S). Used for power-outage outcomes and heavy-feeding events. Fish respiration 400 mg O2/(kg·h) at 26 °C (S: typical range, Timmons & Ebeling).

## 2.6 Light

### Fixture model

Each catalog fixture carries: electrical watts; channels, each with a spectral power distribution (10 nm bins, 380–720 nm) and a PAR map at a reference height (a grid of µmol m⁻² s⁻¹ at 30 cm in air, 10 cm pitch) from manufacturer or independent tests (H: BRS "PAR map" tests, Reef Builders, Sanjay Joshi's measurements). Generic fixtures carry a lambertian-with-optics distribution parameterized by beam angle.

PAR at a point P in water from fixture at height h above the surface:

PAR(P) = Σ_fixtures Σ_channels I_ch(θ) · τ_surface(θ) · exp(−K_d · s_water) / r² · dimmer_ch(t) + PAR_diffuse

where I_ch(θ) is the angular intensity inverted from the reference map, r the distance to the refracted ray's entry point, s_water the in-water path length (Snell's law, n = 1.34), τ_surface the Fresnel transmittance (0.98 at normal incidence to 0.5 at grazing), and K_d the diffuse attenuation. PAR_diffuse = 0.12 × local direct (E; wall reflections and scattering).

| Constant | Value | Tag | Note |
|---|---|---|---|
| K_d | 0.35 m⁻¹ clean water with carbon, 0.6 without carbon (yellowed) | E | Clear ocean water is 0.03–0.05 m⁻¹ at 450 nm (S: Jerlov type I; Morel 1988); tank water carries DOC and particles. Hobby PAR measurements show most of the drop with depth is geometric spread, not absorption (H: BRS tests at 12/18/24 in) |
| Shading | ray-cast against the scape's coarse SDF; fraction of fixture rays reaching the point | | computed per render, cached per colony |

Shown to the user as a PAR overlay in the console.

**As built (catalog 2026.10.4, 2026-10-18).** The first real fixture is the EcoTech Radion XR30 G6 Pro. Per-channel spectra in 10 nm bins and a PAR map in numbers are not published by anyone for it: BRS's XR30 tests (G4 Pro, G5 Blue) give their numbers only in video, and every spectrum found is a graph. So the fixture is built from Sanjay Joshi's independent in-air test (Reef Builders, 1 Oct 2024; LI-COR LI-192SA on a 3 × 3 ft grid at 3" pitch, UPRtek PGN100 spectrometer), in the same form as the generic fixtures, one channel:
- on-axis intensity from his peak readings, 173.2 µmol at 61 cm and 111.7 at 76 cm, which agree on inverse square (716 and 717 at 30 cm);
- the angular term cos^n(θ) fitted to his 2 × 2 ft and 3 × 3 ft averages, n = −0.25 (a wide 126° optic throws slightly more light off-axis; the validator allows −2 to 20);
- band fractions instead of a 10 nm spectrum: 70 % of PAR photons in 400–500 nm (120.7 of 172.0 µmol), the blue fraction the pigment model reads.

Per-channel dimmers wait for per-channel data on the G6 (the only per-channel test found is DaniReef's of the G5 Blue, a different LED mix).

**Validation** (an automated test). The fixture model in air, at 100 %, against Joshi's measurements, within 15 %:

| Joshi 2024, in air | Measured | Model | Difference |
|---|---|---|---|
| Peak, 61 cm | 173.2 | 173.4 | +0.1 % |
| Peak, 76 cm | 111.7 | 111.0 | −0.6 % |
| 2 × 2 ft mean, 61 cm | 150.9 | 148.3 | −1.7 % |
| 3 × 3 ft mean, 61 cm | 124.9 | 130.0 | +4.1 % |
| 2 × 2 ft mean, 76 cm | 103.0 | 99.9 | −3.0 % |
| 3 × 3 ft mean, 76 cm | 92.4 | 90.8 | −1.8 % |

This is an in-air check of the fixture model only; the surface (Fresnel) and water attenuation (K_d, table below) are separate terms with their own sources. An earlier draft checked against BRS at 45 cm; that was dropped because BRS's numbers are not published in a usable form and nobody measured at 45 cm. For reference, the model puts 320 µmol 25 cm deep under an XR30 hung 20 cm above the water (45 cm from the fixture), 171 at 40 cm deep.

### Photoperiod and schedules

Per-channel piecewise-linear schedules, including sunrise and sunset ramps, midday peaks, and "cloud" randomness (user option). Seasonal mode varies photoperiod ±1 h and temperature ±1 °C across a 365-day sim year (S: real reefs; H: a known spawning trigger approach, Craggs et al. 2017 at the Horniman Museum).

### Moon

Lunar phase from the sim date (synodic period 29.530588 d; S). Moonlight intensity for rendering only; its PAR is < 0.01 µmol m⁻² s⁻¹ and contributes nothing to growth (S). Phase is a cue in §2.7 spawning (S: Harrison et al. 1984; Babcock et al. 1986).

### Acclimation and photoinhibition

I_acc relaxes toward the 14-day mean daily PAR: dI_acc/dt = (I_daymean − I_acc)/τ_acc, τ_acc = 14 d (S: photoacclimation over 1–3 weeks, Anthony & Hoegh-Guldberg 2003; Hoogenboom et al. 2006). Photosynthesis uses the Platt et al. 1980 form with inhibition:

P(I) = P_s · (1 − e^(−α I / P_s)) · e^(−β I / P_s) − R

with I_k = P_s/α scaled by acclimation: I_k = I_k,sp · (I_acc / I_ref,sp)^0.5 (E). Photodamage stress rate = k_pi · max(0, I_noon − 1.6 · I_acc − 50) (E; the sudden-jump rule hobbyists follow: new corals under lights at full power bleach, the same lights ramped over weeks do not; H: BRS acclimation guides). UV is not modeled at launch.

## 2.7 Coral biology

### Groups and parameters

| Group | Examples at launch | ext_max cm/yr | G_max mg CaCO3 cm⁻² d⁻¹ | I_k,sp | NO3 window ppm | PO4 window ppm | Tag |
|---|---|---|---|---|---|---|---|
| SPS branching | Acropora, Seriatopora, Pocillopora, Stylophora | 10 (Acropora) 6 (others) | 1.5 | 250 | 2–10 | 0.03–0.10 | S: Dullo 2005 reviews Acropora extension 2–20 cm/yr; Marubini et al. 2001 calcification 0.5–2 mg cm⁻² d⁻¹; windows H |
| SPS plating/encrusting | Montipora capricornis, M. digitata, encrusting Montipora | 8 (radial) | 0.8 | 220 | 2–15 | 0.03–0.15 | S/H as above |
| LPS fleshy | Euphyllia, Trachyphyllia, Fungia, Acanthastrea | 3 | 0.4 | 150 | 5–25 | 0.03–0.25 | H: Borneman; Sprung & Delbeek; S for I_k ranges, Chalker 1981 |
| LPS massive | Favia/Favites, Goniopora, Cyphastrea | 1.5 | 0.3 | 160 | 5–25 | 0.03–0.25 | H |
| Soft | Sarcophyton, Sinularia, Zoanthus, Discosoma, Xenia | 6 (area) | 0 | 120 | 5–40 | 0.05–0.5 | H: Borneman; Sprung & Delbeek |
| Anemones | Entacmaea quadricolor | 8 | 0 | 200 | 2–25 | 0.03–0.25 | H |
| Clams | Tridacna crocea, T. maxima, T. derasa | 4 | 0.5 | 300 | 1–20 | 0.02–0.2 | S: Fitt & Trench 1981; H |

### Growth and calcification

Extension: dL/dt = ext_max/365 · f_I · f_T · f_Ω · f_nut · f_pH · f_allelo · H · (1 − L/L_max)
Calcification: dM/dt = G_max · A(L) · f_I · f_T · f_Ω · f_nut · f_pH · H, with A(L) the species shape's surface area (branching A ≈ 3.2 L²; plating A ≈ 0.8 L²; massive A ≈ 1.6 L²; E: geometric factors)

| Factor | Form | Constants | Tag | Source |
|---|---|---|---|---|
| f_I light | tanh(I_day / I_k) with "light-enhanced calcification" ratio light:dark 3:1 | | S | Chalker 1981 (tanh form); Gattuso et al. 1999 review (LEC) |
| f_T temperature | exp(−((T − T_opt)/σ_T)²), T_opt 26.5 °C, σ_T 3 °C (tank-acclimated; reef studies peak 27–28) | | S/E | Marshall & Clode 2004 (Galaxea calcification peaks near 27–28 °C); Jokiel & Coles 1977 (growth optimum 26–27 °C Hawaii) |
| f_Ω saturation | clamp((Ω − 1)/(Ω_ref − 1), 0, 1.3), Ω_ref = 3.8 | | S | Langdon & Atkinson 2005; Gattuso et al. 1998; Marubini et al. 2008. This is the mechanism behind the hobby's "high alk grows faster" (H) |
| f_nut nutrients | smooth window on NO3 and PO4 (table above); below window: ×0.6 and stress; above: ×(1 − 0.4·excess) | | H / S | Below: Rosset et al. 2017 (phosphate starvation under high N:P destabilizes symbiosis); Wiedenmann et al. 2013. Above: Marubini & Davies 1996 (nitrate reduces calcification); Ferrier-Pagès et al. 2000 (phosphate). **Flag:** laboratory reductions occur at nitrate levels (1–20 µM, 0.06–1.2 ppm) far below hobby "high" values; the window uses hobby consensus and the literature only for direction |
| f_pH | smoothstep(7.6, 7.9, pH_NBS) | | H | Holmes-Farley "Low pH"; mechanistically pH acts through Ω, which is already in f_Ω; this factor only adds the hobby's observed slowdown below 7.8 (E) |
| f_allelo | 1 − 0.3 · SoftIndex · (1 − GAC_active) | | S/H | Coll et al. 1982; Sammarco et al. 1983 (Sinularia/Sarcophyton terpenes reduce scleractinian growth); H: activated carbon mitigates |

Alkalinity demand emerges: ΔTA = −2 · ΔM / 100.09 eq. **Check:** 40 mid-size SPS colonies in a 120 gal with f-factors ≈ 1 consume ≈ 12 g CaCO3/day ≈ 1.5 dKH/day, inside the 1 to 3 dKH/day hobbyists report for SPS-dominated tanks (H: Holmes-Farley "Calcium and alkalinity").

### Health, stress, bleaching

dH/dt = +r_rec · (1 − H) − Σ stress rates, r_rec = 0.02 /d (E; weeks to recover).

| Stressor | Rate | Tag | Source |
|---|---|---|---|
| Thermal | k_T · max(0, T − T_b)², T_b = 30 °C default (species offset ±1); accumulator DD = ∫max(0, T − T_b) dt; bleaching when DD ≥ 1 °C·d (zooxanthellae lost with τ = 2 d), mortality pressure when DD ≥ 6 °C·d (resets at 0.3 °C·d/day below T_b) | E anchored on S | NOAA Coral Reef Watch DHW: bleaching at 4 °C·weeks, mortality at 8 °C·weeks over MMM + 1 °C (Liu et al. 2006). Tanks bleach faster than reefs because O2 drops and there is no thermal mass; the squared term and the shorter accumulation are E |
| Cold | k_C · max(0, 22 − T) | H | Sprung & Delbeek: most reef corals stressed below 22 °C |
| Photodamage | §2.6 | E | |
| Alkalinity level | below 6.5 dKH: 0.03 /d per dKH; above 12 with f_nut low: 0.02 /d per dKH | H | Holmes-Farley recommended range 7–11 dKH; "burnt tips" at high alk with low nutrients is Reef2Reef consensus |
| Alkalinity swing | 0.05 /d per dKH of 24-h range beyond 1 dKH, SPS only | H | Reef2Reef, BRS "alkalinity stability" guidance; no primary source found (a search on 2026-10-18 found none either; Holmes-Farley's alkalinity article tolerates "bouncing" above 7 dKH and gives no rate limit) |
| Salinity | outside 30–38 psu: 0.02 /d per psu; 24-h range > 2 psu: 0.1 /d | H | Sprung & Delbeek |

The 24-h range is a ring of hourly highs and lows over the last 24 sim-hours, updated at the start of every step. Water changes and doses are applied between steps, so a jump they cause counts in full for the next 24 sim-hours; a per-step rate never saw them. Implemented 2026-10-18, after nine back-to-back water changes in a test tank went unpunished: before that, neither swing term was in the code. Measured: nine 10 % changes in a row (10.7 → 9.0 dKH) cost an Acropora about 0.03 health, one 30 % change (0.8 dKH) nothing, a 3 dKH soda-ash dose about 0.10. The console's Water tab previews the swing a change or dose would cause and warns before it is applied.
| Ammonia | f_NH3 · TAN > 0.02 mg/L: 0.5 /d per mg/L | S | Aquaculture chronic thresholds (Timmons & Ebeling) |
| Phosphate starvation | PO4 < 0.02 ppm, whatever the nitrate: 0.01 /d, SPS 0.02 /d; the symbiont target × (0.7 + 0.3 · s(0.005, 0.03, PO4)) | S (direction) / E (sizes) | Rosset et al. 2017 (phosphate deficiency promotes coral bleaching); Wiedenmann et al. 2013 (P limitation under N enrichment); hobby "STN from bottomed-out nutrients". Until 2026-10-18 the stress needed NO3 < 1 ppm as well and the symbionts ignored phosphate, so a carbon-dosed tank at 0.00 PO4 and 1–2 ppm NO3 looked healthy |
| Aggression | contact or sweeper radius: both lose 0.05 /d weighted by hierarchy | H | Borneman "Aquarium Corals"; Sprung & Delbeek: Euphyllia sweepers 5–10 cm, Galaxea 15–30 cm, Favia/Favites night sweepers ~5 cm |
| Pests (optional) | §2.9 | | |

Zooxanthellae and color: dZ/dt = (Z* − Z)/τ_Z, τ_Z = 21 d, with Z* = 1 · g_nut(NO3, PO4) · g_I(I_acc) where g_nut rises to 1.6 at high nutrients and g_I falls from 1 at I_acc = ½·I_ref of the species to 0.6 at 2·I_ref (relative to the species since 2026-10-18; it was 100 to 400 µmol for every species, which paled SPS under SPS light) (S: Dubinsky et al. 1990; Marubini & Davies 1996 for nutrient-driven increase; Titlyanov et al. 2001 for light-driven decrease; the shape is E). Thermal stress drives Z toward 0 while stressed; bleached = Z < 0.35; recovery over 4 to 8 weeks once stress ends (S: typical reef recovery, months; E for the tank constant). Pigment C relaxes toward the blue-light fraction × species pigment capacity with τ = 21 d (S: D'Angelo et al. 2008 blue-light up-regulation of GFP-like proteins). How a colony looks comes from one rule, used both for the console's label and for the coral's colour in the 3D view. What a keeper calls colour is host pigment C, which rises with light, while symbionts are brown and thin out under strong light. With z = Z, v = C / pigment capacity and s a smoothstep: bleached = s(0.55 → 0.25, z); pale = (1 − bleached) · s(0.75 → 0.5, z) (thinned symbionts, washed out); of the rest, the share dull = max(1 − s(0.35 → 0.65, v), s(1.1 → 1.4, z)) is browned (little pigment, or a symbiont bloom from high nutrients) and the remainder colored. The label is the largest; texture weights are colored + 0.4·pale, browned, bleached + 0.6·pale, and skeleton for a dead colony. Thresholds are E: under blue-leaning light (blue fraction 0.6) an Acropora is browned below about 195 µmol acclimated, colored from about 270 and still colored at 470; a hammer or war coral colored from about 135; corals in high nitrate brown out whatever the light. Rendering blends the "colored-up", "browned-out" and "bleached" texture sets by these weights.

### Fragging, mounting, spawning

- Frag: parent L reduced by the frag's share, both get H ← H − 0.2 and a 14-day recovery with growth 0; frag starts at 2–4 cm (H).
- Mounting: a colony encrusts and locks to rock after 30 d; moving it afterward costs H −0.1 (E for the size; recovery follows r_rec, a half-life of about 35 d). Direction and timescale from field transplants (S, fetched 2026-10-18): moved *Acropora palmata* lost more tissue and bleached more at 2 weeks than colonies reattached in place, showed a slight growth and survival penalty at 3 months and none detectable at a year (Forrester et al. 2012, Bull Mar Sci, doi:10.5343/bms.2012.1016); fragments lost a little tissue in their first 3 months (Forrester et al. 2014, doi:10.1002/aqc.2374). Those moves involved breaking and outplanting, so they bound an in-tank re-glue from above; no aquarium measurement was found. The console states the cost before a move. Changing only a colony's shading (exposure) is not a move and costs nothing (until 2026-10-18 it was charged as one).
- Spawning: eligible when age in tank ≥ 365 d sim, L ≥ 0.6 L_max, H > 0.8, T within 26.5–29.5 °C, species flagged as broadcast spawner; on nights 3 to 7 after full moon (S: Babcock et al. 1986: 3–6 nights after the full moon; Harrison et al. 1984), each eligible colony spawns with probability 0.3 per qualifying night if seasonal mode places the date in the species' season, 0.05 otherwise (E). Effects: viewer event (gamete bundles, 20 minutes), DOC pulse of 0.5 mg C/L per colony next day (E; hobby reports of post-spawn water fouling, e.g. Craggs et al. 2017 ex-situ spawning).
- Anemones: wander toward better light/flow when f_I or flow is outside preference for 7 d; a wandering anemone can contact corals (aggression) or a powerhead (death unless guarded; H).

## 2.8 Fish biology

| Process | Model | Constants | Tag | Source |
|---|---|---|---|---|
| Growth | dL/dt = K (L_∞ − L) · f_feed · (0.3 + 0.7 H) | K, L_∞ per species | S where FishBase lists them (e.g. Zebrasoma flavescens K ≈ 0.5/yr, Claisse et al. 2009); otherwise E from max size and reported years to adulthood | |
| Mass | W = a L^b | a, b from FishBase length-weight | S | |
| Ration | 2 % W/day for small planktivores, 1.2 % for tangs, 3 % for juveniles; satiation decays at 0.5 /d | | H/E | Aquaculture rations 1–3 % (S: Timmons & Ebeling) |
| Ammonia | §2.5 | | S | |
| Bioload display | g N/day ÷ V | | | |
| Compatibility | pairwise matrix: conspecific/congeneric aggression (Zebrasoma × Zebrasoma, Centropyge × Centropyge, two clownfish species), predator/prey (lionfish vs small gobies), harassment (sixline wrasse, dottybacks), schooling minimums (Chromis ≥ 5), minimum tank length per species | | H | Scott W. Michael "Reef Aquarium Fishes"; Fenner "The Conscientious Marine Aquarist"; LiveAquaria minimum tank sizes |
| Aggression effect | victim H −0.03 /d per aggressor, halved if tank length ≥ 1.5 × minimum; can kill | | E | |
| Stress | T outside species range: 0.05 /d per °C; pH < 7.7: 0.1 /d; 24-h salinity range > 2 psu (about SG 0.0015; the same window as the corals'): 0.1 /d; f_NH3·TAN > 0.05: 0.5 /d per mg/L; O2 < 4 mg/L: 1 /d; hunger: sat < 0.3: 0.05 /d | | H/S | Ammonia and O2 thresholds S (aquaculture); others H |
| Nitrite | no stress below 1 mg/L N in seawater | | S | Chloride competition reduces nitrite toxicity in marine fish (Lewis & Morris 1986) |
| Nitrate | mild stress above 100 ppm NO3 | | H | Hobby consensus; S for aquaculture tolerance (Timmons & Ebeling) |
| Disease | new fish carry Cryptocaryon with p = 0.3 per purchase (fish bought together share their source system's status; E) reduced to 0.03 by a 30-day quarantine (H); outbreak when stress high; lifecycle 7–14 d at 26 °C (S: Colorni & Burgess 1997); fallow 76 d clears a tank (S: Colorni & Burgess; H: hobby standard); velvet (Amyloodinium) p = 0.1, faster, lethal in 3–5 d untreated (H: Fenner) | | | |
| Cryptocaryon cycle (as built 2026-10-18) | Infected fish shed a free-living stage into the display: load L (0..1), dL/dt = 0.025 · Σ ich · (1 − L) − 0.05 · L; a clean fish in the display is infected once per sim-day with p = 0.3 · L (doubled when stressed); L is cleared once the display has had no fish for 76 days. Hospital tank (`hospital`, all fish or listed ones): copper clears 0.15 of a fish's load a day, costs 0.005 health a day; a 14-day course. `returnFish` puts them back as they are. No infected fish, no load | | E (rates) / S (fallow) / H (course) | Colorni & Burgess 1997 (tomonts develop 3–28+ days at 24–26 °C, up to 72 days cooler: the 76-day fallow); copper course: Fenner, *The Conscientious Marine Aquarist* (hobby standard 14 days at therapeutic level). Before 2026-10-18 ich lived only on fish and spread fish to fish, so removing the infected fish ended an outbreak and the fallow did nothing. Limitation: no acquired immunity, so an untreated infection progresses to death; hardy fish often survive real outbreaks |
| Lifespan | hazard rises after 70 % of max age; max ages from literature where available (Amphiprion > 15 yr; Zebrasoma > 20 yr) | | S/H | Buston & García 2007 (clownfish longevity); H for others |
| Grazing | tangs and blennies (catalog fish `algaeGrazing`) and the cleanup crew (§2.8b) reduce Alg | | H/E | |

## 2.8a Coralline algae (rock surface; as built 2026-10-19)

Crustose coralline algae (CCA) are tracked as one tank-wide coverage of the rock, from 0 (bare dry rock) to 1 (rock crusted over). The coverage drives how the rock looks in the 3D view, together with a seasoned-biofilm index and the nuisance indices. The step runs every integration step, after the nuisances, and is not gated by the pests extension.

dC/dt = µ · max(C, seed) · (1 − C) · f_Ω · f_Mg · f_P · f_I · (1 − Alg) − (d + bleach + o · Alg + diss) · C

- **Seed rule.** Spores arrive three ways: with corals (frags and their plugs carry them), with live rock, and from a bottle (added 2026-10-21). A tank with no corals and C = 0 does not grow coralline, so a fish-only tank on dry rock stays bare. Once C > 0 the seed floor keeps it seeded.
- **Live rock.** It comes in at setup ("Start from scratch on live rock", `newTank({ liveRock: true })`, share 1) or later (`addLiveRock { share }`, the share of the rockwork that is live, default 0.2). Each brings in, scaled by the share, never lowering a value:
  - C to max(C, 0.02 × share, seed), the level at which pale spots show;
  - the nitrifier guilds (aob, aoa, nob) to 0.05 × share, the bottled-bacteria head start;
  - denitCap to 0.1 × share.

  With the pests extension on, it also brings hair algae 0.02 × share and aiptasia with p = 0.3 × share, keyed on the rock's id.
- **Bottled coralline.** `treat coralline_bottle` seeds an unseeded tank at C = seed. Each dose adds 0.05 to a boost b (state `ccaBoost`), capped at 0.15, and growth runs at µ × (1 + b). b halves every 7 days, so a daily dose holds it at the cap.
- **Calcification.** Coverage does not draw on Ca or alkalinity yet. Coralline calcification is still the fixed term in §2.4 (`coralDemandGCaCO3PerDay`).
- **Migration.** A tank saved before this field existed loads crusted at C = 0.85 if it has corals and denitCap > 0.5 (matured rock, as the running-reef preset starts). Otherwise it loads at 0. The running-reef preset starts at 0.85.
- **Equilibrium.** In the default tank (dKH 8, Ca 420, 550 ppm room, Ω_arag ≈ 3.9, PAR at the rock ≈ 167 µmol m⁻² s⁻¹), pale spots (C ≈ 0.02) show at about 7 weeks, colour (C ≈ 0.15) at about 4½ months, half the rock at about 7 months, and coverage settles at about 0.92. At common hobby minimum targets (dKH 9, Ca 400, Mg 1,300) the times are about 6 weeks, 4 months and 6 months; at Ca 440 they are about 5½ weeks, 3½ months and 5½ months.
- **Targets (2026-10-20).** Tuned to an experienced keeper's working targets: calcium above 400 (often 440), dKH 9+ and Mg 1,300+, and a slow grower. Alkalinity and calcium act through Ω_arag rather than as separate ramps, so pH (the house's CO2) counts too and nothing is counted twice.

| Term | Model | Constants | Tag | Source |
|---|---|---|---|---|
| Growth | logistic, µ · C · (1 − C) | µ = 0.035 /d, seed = 0.006 | H | H (hobby timelines from arcreef.com and forum threads): pale spots at 4–8 weeks, pink to purple by 2–4 months, heavy cover at 6–12 months. This is the right order for marginal extension of 0.9–2.3 mm/month (S: Adey & Vassar 1974, St Croix) and 2.9–3.9 mm/month (S: Matsuda 1989, Ishigaki) from spores a few cm apart |
| Saturation | f_Ω = min(1.3, max(0, (Ω_arag − 1.2) / (4.25 − 1.2))); dissolution 0.02 /d per unit Ω below 1.5 | | S (direction and shape), E (low end, cap), H (f_Ω = 1 point) | Kuffner et al. 2008, Nat Geosci 1:114: CCA recruitment and cover fall steeply at elevated CO2. Full end: Ω_arag of dKH 9 and Ca 400 in a 600 ppm house. dKH 7 / Ca 380 grows at ≈ 50 %, and in a 1,000 ppm house at ≈ 20 %. Above that point growth keeps rising with Ω, so more calcium means faster growth (hobby experience; S for the shape: Comeau et al. 2013, Limnol Oceanogr, eight reef calcifiers including corallines; citation still to be verified, calcification rising with Ω well above saturation): Ca 440 runs at ≈ 1.14. It is capped at 1.3 (Ω ≈ 5.2), just below where abiotic precipitation starts (Ω_c 5.5, §2.4) |
| Magnesium | f_Mg = ramp(Mg, 1100, 1300 mg/L) + 0.2 · ramp(Mg, 1300, 1450) | | H | Holmes-Farley: crusts are high-Mg calcite, and low Mg upsets the Ca/alk balance. 1 at 1,300+; then a modest boost to 1.2 at 1,450 and flat beyond, as reefers raise Mg to 1,400–1,500 to push coralline |
| Phosphate | f_P = 1 / (1 + (PO4 / K)²) | K = 0.15 ppm PO4 (0.049 mg P/L; the code converts `po4MgL` × 3.07) | H | Hobby rule of thumb: slow above 0.1 ppm, stalled by 0.3. No experimental source confirmed (§2.13) |
| Light | f_I = I / (25 + I), with I the daily-mean PAR at 60 % of water depth; bleaching loss 0.03 /d per 350 µmol above 350 | | E; S (qualitative) | Burdett et al. 2014, BMC Plant Biol 14:139: red corallines photoinhibit and may bleach at high irradiance |
| Competition | growth × (1 − Alg); overgrowth loss 0.02 /d at Alg = 1 | | S (direction), H | Kuffner et al. 2008: fleshy algae rose as crusts fell |
| Background loss | d = 0.002 /d (grazing, scraping) | | E | Sets mature coverage |
| Seasoned biofilm (derived, no state) | B = min(1, max(aob, aoa, nob, denitCap) / 0.25) | | E | The brown film grows with the bacteria |
| Live rock | C 0.02, nitrifiers 0.05, denitCap 0.1, algae 0.02, aiptasia p 0.3, each × share | | H (coralline, nitrifiers, aiptasia), E (denitCap, algae) | H: cured live rock arrives with scattered coralline patches, shortens the cycle (the reason reefers buy it), and is the classic aiptasia source |
| Bottled coralline | growth × (1 + b); b += 0.05 per dose, cap 0.15, half-life 7 d; seeds at C = seed | | E | Kept small on purpose: the hobby evidence that bottled spores speed growth is anecdotal, with no controlled source |

## 2.8b Cleanup crew and display copper (as built 2026-10-23)

Invertebrates are kept by headcount, one group per species. Species values are in the catalog (`INVERTS`, catalog 2026.10.6, per-field provenance). Each group has a condition c in 0..1. The species set and the values are provisional defaults.

dc/dt = 0.1 · fed · (1 − c) − (1 − fed) / starveDays − stress

- **Diet and fed.** An animal is fed by the best of its diets (omnivores switch). Each food is shared among everyone who eats it.
  - Grazing: supply V · (biofilm / 7.6 L + algae / 1.5 L) in "one Trochus" units, against the sum of grazeUnits × count. Biofilm is the seasoned-rock index of §2.8a. A mature tank's film alone feeds about one snail per 2 US gal.
  - Scavenging: supply is the uneaten food plus 5 % of what the fish ate (waste and detritus), against scavengeGPerDay × count.
  - Aiptasia: fed in full at an aiptasia index of 0.02 and above.
- **Stress** (per day):
  - 0.05 per °C outside the species' range;
  - 0.2 outside its salinity range;
  - 0.15 while the 24-hour salinity swing exceeds its limit;
  - 0.02 for moulting crustaceans below 1,150 mg/L Mg;
  - NH3 and O2 as for fish;
  - copper: 0.1 above a fifth of the species' lethal level, and 1.0 × Cu / lethal at or above it.
- **Deaths.** One animal at a time, at a hazard of 2 × (0.5 − c) × 2 per animal per day below c = 0.5, plus 1/30 a day past the species' lifespan. All die at c = 0. A hermit crab kills a snail at snailPredation × (2 − fed) per crab per day, because the model has no spare shells.
- **Effects on the tank.**
  - Grazing removes 0.004 of the algae index per unit per day in 450 L (scaled by volume), added to the fish's grazing (§2.9). Thirty Trochus graze like one adult yellow tang.
  - Scavenged food is eaten before it rots. It goes to ammonia like the fish's food, not to dissolved organic carbon.
  - Aiptasia eaters add aiptasiaPerDay × count × c to the aiptasia loss.
  - Urchins scrape coralline at corallineGrazing × count × c × 450/V, added to the coralline loss (§2.8a).
  - Cleaner shrimp add parasitePick × count × c (capped at 0.1/day) to the loss of the free-living ich load. They eat tomonts and barely touch the feeding stage on fish (Vaughan et al. 2018), so they slow an outbreak but do not cure it.
- **Treatments that became animals.** `treat peppermint` now adds two peppermint shrimp, and `treat berghia` adds five Berghia. A tank treated before 2026-10-23 keeps its older aiptasia-control setting.
- **Not replaced.** The fixed coralline background loss (grazing snails and urchins) and the fish grazing stay as they are, so a tank without crew runs exactly as before. The crew adds on top.
- **Display copper.** `treat copper_display` holds dissolved copper at 0.2 mg/L for 14 days, the hobby's ionic copper course.
  - At the therapeutic level it clears ich from the display fish at 0.15/day (as the hospital tank) and takes the free-living load down an extra 0.3/day.
  - It costs fish 0.005/day and corals 0.5 × Cu per day, so corals are dead in about 10 days.
  - After the course it falls 5 % a day, plus 5 % a day with carbon, and water changes dilute it. That leaves too much for invertebrates (over 0.01 mg/L, Yanong FA165) for 3–8 weeks.
  - This is the trade-off against the hospital tank (§2.8), which keeps copper away from the display.

| Constant | Value | Tag | Source |
|---|---|---|---|
| Film grazing supply | 1 unit per 7.6 L at biofilm 1 | H/E | Hobby stocking (1 snail per 1–2 gal assumes some algae; overstocked snails starve) |
| Algae grazing supply | 1 unit per 1.5 L at algae 1 | E | |
| Grazing per unit | 0.004 /d of the algae index in 450 L | E | Sized against the yellow tang's 0.12 /d (catalog) |
| Detritus share | 5 % of the food eaten | E | |
| Aiptasia fed at | index 0.02 | E | |
| Cleaner cap | 0.1 /d on the free-living load | E | Direction: Militz & Hutson 2015 (S), Vaughan et al. 2018 (to verify) |
| Recovery, temperature, salinity, NH3, O2 stress | as the fish model | E | |
| Swing stress | 0.15 /d | H/E | Invertebrates are drip-acclimated because swings kill them |
| Moult stress | 0.02 /d below 1,150 mg/L Mg | E | Hobby lore; no source found |
| Copper course | 0.2 mg/L for 14 d | H | UF/IFAS FA165: 0.15–0.20 mg/L free copper |
| Invertebrates safe at | ≤ 0.01 mg/L | S | Yanong, UF/IFAS FA165 |
| Copper loss after the course | 5 %/d, +5 %/d with carbon | E | Residual copper in rock is a known hobby problem |
| Coral damage | 0.5 × Cu /d | E | H: copper kills corals |

## 2.9 Nuisances and events (optional extensions, off by default)

| Extension | Trigger model | Resolution the user can apply | Tag |
|---|---|---|---|
| Green hair algae | dAlg/dt = µ_a · Alg · f_I · NO3/(1+NO3) · PO4/(0.03+PO4) − grazing − manual removal | grazers, nutrient export, GFO, lights | H |
| Cyanobacteria | grows where flow < threshold and DOC high and PO4 > 0.05; bloom event | flow, feeding, water changes, chemiclean (modeled as antibiotic: also sets B_AOB −30 %) | H |
| Dinoflagellates | probability rises when NO3 < 1 and PO4 < 0.02 for > 14 d (post-GFO or carbon-dosing crash); growth under light | raise nutrients, UV, blackout, dosing bacteria | H (Reef2Reef consensus; no primary aquarium source found) |
| Aiptasia | arrives with frags p = 0.15; grows with feeding | peppermint shrimp, Berghia, spot treatment | H |
| Bryopsis, vermetids, flatworms | later | | |

Implementation: each nuisance is an index in [0, 1] with logistic growth, a spore floor while the extension is on, and the drivers above as saturating factors: hair algae µ = 0.4/day at saturating light (half-saturation 80 µmol at the rock surface) with uptake of N and P at algal tissue ratios (S: Atkinson & Smith 1983) and grazing per tang or blenny from each fish's catalog `algaeGrazing` (E), scaled by the fish's size; cyanobacteria µ = 0.6/day gated by turnover below 15×/h, DOC and PO4 above 0.05 ppm, fading at 0.08/day; dinoflagellates have an outbreak hazard that rises after 14 days of NO3 < 1 and PO4 < 0.02 (2 %/day per further fortnight, capped at 10 %), then grow at 0.5/day under light, lose 0.12/day with competition restored, and 0.25/day per tank-turnover-per-hour of UV; aiptasia arrives with p = 0.15 per frag (keyed on the frag id so replays agree) and grows at 0.1/day under heavy feeding. Effects: algae smother, cyano cover, dino toxicity and aiptasia stings cost coral health at 0.01, 0.01, 0.03 and 0.01 per day at index 1; dinos stress fish at 0.02/day. Treatments: scrub (−70 % algae), siphon (−70 % cyano), chemiclean (cyano to 1 %, nitrifiers −30 %), blackout (lights off for N days, affects corals too), aiptasia spot (−50 %), peppermint shrimp (0.1/day), Berghia (0.25/day). Equipment failures draw per step against the MTBF table with the equipment-failure RNG stream: heaters stick on or off (50/50), pumps stop (`failed`, cleared by re-applying the pump list), the ATO float sticks and adds 8 % of the volume as fresh water, a dosing pump delivers two days' dose at once, a fixture goes dark, and the power fails for 1–12 h at 0.5/yr and comes back by itself. All these numbers are E or H, and open to challenge.
| Equipment failure | MTBF: heater 3 yr (stuck on or off, 50/50), powerhead 5 yr, return 6 yr, ATO sensor 2 yr (flood → salinity drop), dosing pump runaway 4 yr, light 6 yr | replace | E |
| Power outage | user-triggered or random (0.5/yr, 1–12 h): T relaxes to room, O2 declines at respiration rate, pH drifts to air equilibrium; fish die when O2 < 2 mg/L | battery backup air pump, generator | S for physics; outage frequency E |
| RODI TDS creep | top-off adds PO4 at user-set ppm | change membranes/resin | H |
| Coral spawning | §2.7 | | |

## 2.10 Equipment models

| Equipment | Inputs from catalog | What it does in the model |
|---|---|---|
| Light fixture | watts, channels with SPD, PAR map at 30 cm, beam geometry | §2.6; heat η_L · watts · dimmer |
| Heater | watts | thermostat, hysteresis ±0.3 °C, failure modes |
| Chiller | BTU/h, COP | capacity when above setpoint |
| Return pump | watts (submerged or external), flow L/h | heat, surface agitation (k_surface), sump exchange |
| Wavemaker/powerhead | watts, flow L/h, pattern | heat, flow field for cyano/detritus/coral flow preference, k_surface |
| Skimmer | rated volume, air draw L/h, watts | DOC export fraction, k_skimmer (CO2), heat |
| Dosing pump | channels, mL/day each, solution | §2.4 per solution |
| Kalk reactor/stirrer | in ATO line, saturation 0–100 % | §2.4 |
| Calcium reactor | effluent mL/min, effluent dKH, media type | §2.4 |
| ATO | on/off, source TDS | holds V; adds P at TDS creep |
| Refugium | volume, light PAR, macroalgae mass | §2.5 uptake; reverse photoperiod stabilizes pH (emergent from DIC) |
| Media reactor GFO / GAC | grams, change interval | PO4 binding; K_d lowering; allelopathy removal |
| UV sterilizer | watts, flow | reduces free-swimming parasite stage and algae spores (E) |
| Fans / lid | fan on/off; lid | v and E in evaporation |
| Controller (Apex/GHL) | later | reads real tank into the model ("mirror") |

## 2.11 Species parameter schema

Catalog entries are JSON validated against a schema. Every numeric field is an object `{ "v": value, "tag": "S|H|E|T", "src": "ref-id", "note": "" }`. The 3D renderer reads the same file.

```json
{
  "id": "acropora_millepora",
  "kind": "coral",
  "taxonomy": { "genus": "Acropora", "species": "millepora", "common": "Millepora" },
  "group": "sps_branching",
  "morphology": { "form": "branching_corymbose", "areaFactor": {"v": 3.2, "tag": "E"} },
  "size": { "fragCm": {"v": 3}, "maxCm": {"v": 30, "tag": "H", "src": "borneman2001"} },
  "growth": { "extMaxCmPerYr": {"v": 10, "tag": "S", "src": "dullo2005"}, "gMax": {"v": 1.5, "tag": "S", "src": "marubini2001"} },
  "light": { "ik": {"v": 250, "tag": "S", "src": "chalker1981"}, "iRef": {"v": 300, "tag": "E"}, "inhibit": {"v": 700, "tag": "E"} },
  "temperature": { "optC": {"v": 27}, "sigmaC": {"v": 2.5}, "bleachC": {"v": 30, "tag": "E"} },
  "nutrients": { "no3Ppm": [2, 10], "po4Ppm": [0.03, 0.10], "tag": "H", "src": "r2r-sps-nutrients" },
  "aggression": { "sweeperCm": {"v": 0}, "rank": 3, "allelopathic": false },
  "reproduction": { "broadcast": true, "seasonMonths": [10, 11], "src": "babcock1986" },
  "pigments": { "capacity": 0.9, "blueResponse": 0.8, "textureSets": ["colored", "browned", "bleached", "skeleton"] },
  "assets": { "model": "corals/acropora_millepora", "lods": 3, "growthStages": [3, 8, 15, 25] },
  "catalogVersion": "2026.10"
}
```

Fish entries carry `vonBertalanffy {k, lInf}`, `lengthWeight {a, b}`, `rationPct`, `dietProtein`, `ranges {tempC, pH, sg}`, `minTankCm`, `compat` rules, `behavior`, `lifespanYr`, `diseaseSusceptibility`, and `assets {model, rig, animations}`.

Equipment entries carry the fields in §2.10 plus `brand`, `model`, `msrpUsd`, `dataSource`.

## 2.12 Validation scenarios

Each scenario is a fixture file (initial state + schedule), run headless by the test suite on every catalog or model change. Tolerances are generous where the hobby evidence is qualitative.

| # | Scenario | Setup | Expected outcome | Tolerance |
|---|---|---|---|---|
| 1 | New-tank cycle | 120 gal, 120 lb dry rock, 2 ppm TAN dosed day 0, no bacteria added | TAN peak then NO2 peak 7–14 d later; both < 0.1 by day 18–35; NO3 ends 5–10 ppm | days ±5 |
| 2 | Bottled bacteria | as 1 with B₀ = 0.1 | cycle complete < 10 d | ±3 d |
| 3 | Fish added to uncycled tank | 120 gal, 3 fish, feeding, no cycle | TAN > 0.5 by day 4; un-ionized NH3 > 0.05; fish health falls; first death by day 10–20 | implemented, passing (adult yellow tang + two clowns, 5 g/day: TAN 1.1 and NH3 0.11 on day 4, first death day 11.4; the same fish in a cycled tank stay healthy. Fish share one stress term, so all three die together, and survivors recover once the nitrifiers catch up on day 12–13) |
| 4 | No top-off | 120 gal open top, fan on, ATO off | SG rises ≈ 0.0005/day; corals stressed after SG > 1.030 (≈ day 10–14) | ±30 % on rate; implemented, passing (from SG 1.025: ≈ 0.0006/day over the first 8 days, speeding up as the volume shrinks; SG passes 1.030 on day 9; Acropora health 1.0 until SG 1.028 (coral salinity stress starts at 38 psu), then 0.93, 0.78 and 0.40 on days 9, 11 and 14 against 1.0 with the ATO on. From the default 35 psu start SG passes 1.030 on day 6) |
| 5 | House CO2 | 8 dKH, skimmer on, room 420 vs 1,000 ppm | pH 8.24 vs 7.93 ± 0.07 NBS; outdoor air line raises pH ≥ 0.15 within 24 h | implemented, passing |
| 6 | Kalk vs two-part vs reactor | identical SPS tanks, each held at 8 dKH by its method | kalk tank pH +0.15 to +0.25 above two-part; reactor tank −0.1 to −0.2 below | ±0.05 |
| 7 | Snowstorm | kalk overdose to 14 dKH, Ca 520, pH 8.6 | abiotic precipitation lowers alk and Ca within hours; cloudiness event; Ω > 8 | implemented, passing (Ω_ar ≈ 19; ≈ 3 dKH/day at onset, self-limiting within a day) |
| 8 | Low magnesium | Mg 950, alk 9, Ca 450 | alk and Ca fall faster than demand explains; precipitation event; raising Mg stops it | implemented, passing (≈ 0.1 dKH/day with no demand; the same water at Mg 1300 holds; dosing Mg to 1300 cuts the rate > 10×) |
| 9 | SPS alk demand | 40 established SPS colonies in 120 gal, no dosing | alk falls 1–2 dKH/day; growth stalls below 6.5 dKH; health falls below 6 |implemented, passing |
| 10 | Lights too strong | fresh Acropora frags placed under 500 µmol with I_acc 150 | bleaching within 7–14 d; same frags ramped over 4 weeks do not bleach |implemented, passing |
| 11 | SPS in high nutrients | NO3 40, PO4 0.3 | Z → browned within 3 weeks; growth 40–60 % of clean reference; no mortality |implemented, passing |
| 12 | Ultra-low nutrients | GFO + carbon dosing to NO3 0.2, PO4 0.01 | pale corals, SPS STN events within 3–6 weeks; dino risk event |implemented, passing |
| 13 | Heat wave | 31 °C for 48 h, then 26 °C | bleaching started; most colonies recover over 4–8 weeks; sensitive SPS mortality 10–30 % |implemented, passing (tank-scale thresholds: bleaching at 1 °C·d, zooxanthellae loss τ 2 d, E) |
| 14 | Power outage | 8 h outage, heavy fish load, winter room 18 °C | O2 below 4 mg/L by hour 4–6, fish losses; light load survives |implemented, passing |
| 15 | Spawning | mature Acropora, seasonal mode, full moon | spawning on nights 3–7 after full moon with p ≈ 0.3 per night; none when immature |implemented, passing |
| 16 | Quarantine | add 5 fish with vs without 30-day QT, 100 runs | outbreak incidence ≈ 30 % vs ≈ 3 % | ±10 points; implemented, passing (31 % vs 3 % over the same 100 seeded purchases; outbreak = any fish at ich ≥ 0.5, or lost, 7 days after adding. The p = 0.3 / 0.03 of §2.8 now applies per purchase: fish added in the same step share one carrier draw (E). Read per fish, as before, it gave 80 % vs 16 %) |
| 17 | Heater sizing | 180 gal, 100 W heater, room 18 °C, winter | tank cannot hold 25.5 °C (heater 100 % duty, T falls); 500 W holds |implemented, passing |
| 18 | Determinism | advance 30 d in 1 call vs 720 calls | bit-identical state | exact |

Scenarios 1, 5, 6, 9 and 18 are the core credibility checks and are held to tight tolerances; the rest can be qualitative.

## 2.13 Claims that could not be fully sourced

These are listed here so they can be challenged. If you have a source for any of them, please open an issue.

- The effective nitrifying surface area of live rock (0.5 m²/kg) and sand.
- Coralline algae: the phosphate half-saturation (0.15 ppm PO4), the bleaching light threshold, and the saturation ramp's low end, and the full-growth end points for Ω_arag and Mg (§2.8a, hobby values). Live rock's seeding values and the bottled-coralline boost are hobby values or estimates; the bottle has no source at all.
- Cleanup crew (§2.8b): no measured grazing, scavenging or aiptasia-eating rate was found for any of the seven animals. Those rates are E, sized against the model's existing grazers and the old treatment rates. The copper lethal levels per taxon are E, ranked from secondary summaries. Peppermint shrimp (Rhyne et al. 2004) and Berghia (Monteiro et al. 2020) rates are in paywalled papers and are the first to upgrade.
- Denitrification capacity per kg of rock and the 6-month maturation constant.
- Live rock and sand phosphate sorption capacity and rates.
- The skimmer CO2 equilibration efficiency ε and the surface gas-transfer mapping to flow settings.
- The thermal stress exponent and tank-scale accumulation thresholds (NOAA DHW is a reef-scale index).
- The nutrient optimum windows per coral group (hobby consensus; laboratory effects appear at much lower concentrations).
- Alkalinity swing stress (widely believed, no primary study found).
- Dinoflagellate outbreak triggers (Reef2Reef consensus only).
- GFO capacity per gram, kalkwasser phosphate precipitation yield.
- Light heat fraction into water per fixture type.
- Fish disease introduction probabilities with and without quarantine.
- Equipment MTBFs.
- Spawning probability per qualifying night in captivity.

## 2.14 References

Anthony KRN, Hoegh-Guldberg O (2003) Functional Ecology 17:246–259. Forrester GE et al. (2012) Bull. Mar. Sci., doi:10.5343/bms.2012.1016. Forrester GE et al. (2014) Aquatic Conserv., doi:10.1002/aqc.2374. Bagchi S et al. (2014) PLoS ONE, doi:10.1371/journal.pone.0113515. Babcock RC et al. (1986) Marine Biology 90:379–394. Berner RA (1975) Geochim. Cosmochim. Acta 39:489–504. Borneman E (2001) *Aquarium Corals*, TFH. Buston PM, García MB (2007) J. Fish Biology 70:1710–1719. Chalker BE (1981) Marine Biology 63:135–141. Chave KE (1954) J. Geology 62:266–283. Claisse JT et al. (2009) Marine Ecology Progress Series 389:245–255. Coll JC et al. (1982) Marine Biology 68:XXX (soft coral allelopathy; see Sammarco 1983). Colorni A, Burgess P (1997) Aquarium Sciences and Conservation 1:217–238. Craggs J et al. (2017) Ecology and Evolution 7:11066–11078. D'Angelo C et al. (2008) Marine Ecology Progress Series 364:97–106. Dickson AG (1990) Deep-Sea Research 37:755–766. Dickson AG, Sabine CL, Christian JR (2007) *Guide to best practices for ocean CO2 measurements*, PICES Special Publication 3. Dubinsky Z et al. (1990) Proc. R. Soc. B 239:231–246. Dullo WC (2005) Facies 51:33–48. Emerson K et al. (1975) J. Fish. Res. Board Can. 32:2379–2383. Feldman KS et al. (2010, 2011) Advanced Aquarist, "Total organic carbon" and skimmer performance series. Fenner R (2008) *The Conscientious Marine Aquarist*, 2nd ed. Ferrier-Pagès C et al. (2000) Coral Reefs 19:103–113. Fitt WK, Trench RK (1981) Biological Bulletin 161:213–235. Garcia HE, Gordon LI (1992) Limnol. Oceanogr. 37:1307–1312. Gattuso JP et al. (1998) Global and Planetary Change 18:37–46. Gattuso JP, Allemand D, Frankignoulle M (1999) American Zoologist 39:160–183. Harrison PL et al. (1984) Science 223:1186–1189. Holmes-Farley R (2007) Reefkeeping Magazine 2007-02 (ammonia), reefkeeping.com/issues/2007-02/rhf. Holmes-Farley R, Reefkeeping Magazine and Advanced Aquarist series: "What is alkalinity?", "The relationship between alkalinity and pH", "Low pH: causes and cures", "Calcium and alkalinity", "Magnesium in reef aquaria", "Solving calcium and alkalinity problems", "What your grandmother never told you about lime", "Calcium reactors", "Phosphate: what is it and why should you care", "Iron oxide hydroxide (GFO) phosphate binders", "Specific gravity: oh how complicated". Hoogenboom MO et al. (2006) Coral Reefs 25:xxx. Hovanec TA, DeLong EF (1996) Appl. Environ. Microbiol. 62:2888–2896. Hovanec TA et al. (1998) Appl. Environ. Microbiol. 64:258–264. Jung MY et al. (2022) ISME J, doi:10.1038/s41396-021-01064-z. Incropera FP, DeWitt DP, *Fundamentals of Heat and Mass Transfer*. Jokiel PL, Coles SL (1977) Marine Biology 43:201–208. Langdon C, Atkinson MJ (2005) J. Geophys. Res. 110:C09S07. Lewis WM, Morris DP (1986) Trans. Am. Fish. Soc. 115:183–195. Liu G et al. (2006) Eos 87:xxx (NOAA Coral Reef Watch DHW). Lueker TJ, Dickson AG, Keeling CD (2000) Marine Chemistry 70:105–119. Martens-Habbena W et al. (2009) Nature 461:976–979, doi:10.1038/nature08465. Marshall AT, Clode P (2004) Coral Reefs 23:218–224. Marubini F, Davies PS (1996) Marine Biology 127:319–328. Marubini F et al. (2001) Marine Ecology Progress Series 220:153–162. Marubini F et al. (2008) Coral Reefs 27:491–499. Michael SW (2005) *Reef Aquarium Fishes*, TFH. Millero FJ (1995) Geochim. Cosmochim. Acta 59:661–677. Millero FJ et al. (2001) Marine Chemistry 76:xxx (phosphate adsorption on calcite). Morel A (1988) J. Geophys. Res. 93:10749–10768. Mucci A (1983) American Journal of Science 283:780–799. Neori A et al. (2004) Aquaculture 231:361–391. Platt T, Gallegos CL, Harrison WG (1980) J. Marine Research 38:687–701. Rasheed M et al. (2002) Coral Reefs 21:xxx. Salem S et al. (2006) Biotechnol. Bioeng., doi:10.1002/bit.20822. Sauder LA et al. (2011) PLoS ONE, doi:10.1371/journal.pone.0023281. Redfield AC (1958) American Scientist 46:205–221. Risk MJ, Muller HR (1983) Bull. Marine Science 33:xxx. Rosset S et al. (2017) Marine Pollution Bulletin 118:180–187. Sammarco PW, Coll JC, La Barre S, Willis B (1983) Coral Reefs 1:173–178. Sharqawy MH, Lienhard JH, Zubair SM (2010) Desalination and Water Treatment 16:354–380. Sprung J, Delbeek JC (1994–2005) *The Reef Aquarium* vols 1–3, Ricordea/Two Little Fishies. Tchobanoglous G et al. (2014) *Wastewater Engineering: Treatment and Resource Recovery* (Metcalf & Eddy), 5th ed. Timmons MB, Ebeling JM (2010) *Recirculating Aquaculture*, 2nd ed. Titlyanov EA et al. (2001) Marine Ecology Progress Series 222:xxx. Uppström LR (1974) Deep-Sea Research 21:161–162. Wanninkhof R (1992) J. Geophys. Res. 97:7373–7382. Weiss RF (1974) Marine Chemistry 2:203–215. Wiedenmann J et al. (2013) Nature Climate Change 3:160–164. Zeebe RE, Wolf-Gladrow D (2001) *CO2 in Seawater: Equilibrium, Kinetics, Isotopes*, Elsevier. BRS (Bulk Reef Supply) PAR map tests, heater guide, alkalinity stability guides; Reef Builders fixture tests; Sanjay Joshi lighting measurements (manjoshi.com). Page numbers marked xxx have not been verified yet.
