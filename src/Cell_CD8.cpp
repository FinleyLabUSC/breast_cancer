#include "../inc/Cell.h"
#include "../inc/ModelUtil.h"
#include <cmath>

void Cell::initialize_CD8_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
    state = 6;

    mu = cellParams[0][2];
    kc = cellParams[1][2];
    damping = cellParams[2][2];
    maxOverlap = cellParams[3][2]*cellParams[4][2];
    radius = cellParams[4][2]/2.0;
    deathProb = cellParams[5][2];
    migrationSpeed = cellParams[6][2];
    next_migrationSpeed = migrationSpeed;

    baseKillProb = cellParams[7][2];
    infScale = cellParams[8][2];
    influenceRadius = cellParams[9][2];
    migrationBias = cellParams[10][2];
    divProb_base = cellParams[11][2];
    divProb = divProb_base;
    deathScale = cellParams[12][2];
    migScale = cellParams[13][2];

    migration_speed_base =migrationSpeed;
    kill_prob_base = baseKillProb;
    killProb = baseKillProb;

    location_history.push_back(x);

    death_prob_base = deathProb;

    rmax = 1.5*radius*2;

    init_time = init_tstamp;
}

/**
 * This function simulations inhibition or suppression of CD8+ T cells via PD1-PDL1 binding.
 * First, the minimum between the available PD1 and PDL1 are determined. This assumes a 1-1 binding occurs,
 * thus if there are x PD1 molecules and y PDL1 molecules, with x < y, only x binding "events" will occur.
 * Then, the fraction of PD1 that is bound is calculated. This strength of the inhibition via the binding is directly
 * proportional to this fraction, and occurs with a probability equal to this fraction.
 *
 * @param otherX location of the cell in interacting cell (in direct contact)
 * @param otherRadius radius of the interacting cell
 * @param otherpdl1 level of PD-L1 expressed by the interacting cell
 * @param dt
 * @param master_rng
 * @param temporary_rng
 */
void Cell::cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    // inhibition via direct contact, simulates the binding of PD-1 and PDL1.
    if(state != 6){return;}

    double distance = calcDistance(otherX);
    if(distance <= radius+otherRadius){
        // the line below assumes a one to one binding of pd1 and pdl1. You can only bind the minimum num of molecules of either pd1 or pdl1.
        double percent_binding = std::min(otherpdl1, pd1_available) / pd1_expression_level; // this turns it into a fraction.
        double rnd = master_rng.uniform(0,1,temporary_rng);
        if(rnd < percent_binding){
            // The effect of the pd1-pdl1 binding is a function of the percent_binding and the inhibitory_effect_of_binding.
            next_killProb = next_killProb * percent_binding * inhibitory_effect_of_binding_PD1_PDL1;
            next_migrationSpeed = next_migrationSpeed * percent_binding * inhibitory_effect_of_binding_PD1_PDL1;
            next_death_prob = next_death_prob * (1 + percent_binding) * inhibitory_effect_of_binding_PD1_PDL1;
        }
    }
}

void Cell::cd8_update_properties_indirect() {
    if(state != 6){return;}
    // this assumes that cancer cells don't secrete any factors that affect the CD8's in a contact independent manner
    // posInfluence is M1 + Th
    // negInfluence is M2 + Treg + MDSC
    double posInfluence = 1 - (1 - influences[1])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[5])*(1 - influences[10]);
    double scale = posInfluence - negInfluence;

    // Adjust the properties according to the net signal.
    next_killProb = next_killProb*pow(infScale, scale);
    next_migrationSpeed = next_migrationSpeed*pow(migScale, scale);
}


double Cell::cd8_setProliferationScale(double anti_ctla4_concentration) {
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[8]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[10]);

    // If the cell is severely exhausted, then the drug has no effect
    double anti_ctla4_effect = Hill_function(anti_ctla4_concentration,anti_CTLA4_IC50,anti_CTLA4_hill_coeff) * sensitivity_to_antiCTLA4();

    // Activated CD8 T cells upregulate CTLA-4, and therefor may be more sensitive to its effects. If the cell is more "memory" like then exert a greater effect.
    double scale_anti_CTLA4_effect = (divProb > divProb_base) ? 1.1 : 1;

    // this prevents unexpected negatives in the scale calculation.
    double effective_antiCTLA4_effect = (anti_ctla4_effect * scale_anti_CTLA4_effect < 1) ? anti_ctla4_effect * scale_anti_CTLA4_effect : 1;
    double scale = posInfluence * (1 - influences[5]* (1-effective_antiCTLA4_effect)) - negInfluence;
    return scale;
}


void Cell::cd8_pd1_expression_level(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug) {
    if (type == 3) {
        double temp;
        double influence = 1 - (1-influences[2]) * (1-influences[3])*(1-influences[5])*(1-influences[10]); // interactions with m2, cancer, treg, mdsc
        if (influence >= threshold_for_pd1_induction ) {
            double pd1_increase_amount = (influence - threshold_for_pd1_induction) * pd1_induction_rate * dt;
            pd1_expression_level = (pd1_expression_level < max_pd1_level)?pd1_expression_level : max_pd1_level;
        } else {
            double pd1_decrease_amount = pd1_decay_rate * dt;
            pd1_expression_level -= pd1_decrease_amount;
            pd1_expression_level = (pd1_expression_level > 0) ? pd1_expression_level : 0;
        }

        // determines if the cell is "sensitive" to the drug based on how exhausted it is.
        double fraction_pd1_bound_by_drug = sensitivity_to_antiPD1(anti_pd1_concentration,binding_rate_pd1_drug);

        pd1_drug_bound = pd1_expression_level * fraction_pd1_bound_by_drug;
        pd1_available = pd1_expression_level * (1-fraction_pd1_bound_by_drug);
    }
}


double Cell::Hill_function(double concentration, double EC50, double n) {
    double effect = 1.0 / (1.0 + pow(EC50 / concentration,n));
    return effect;
}

