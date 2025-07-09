#include "../inc/Cell.h"

void Cell::initialize_NK_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
    state = 8;

    mu = cellParams[0][4];
    kc = cellParams[1][4];
    damping = cellParams[2][4];
    maxOverlap = cellParams[3][4]*cellParams[4][4];
    radius = cellParams[4][4]/2.0;
    deathProb = cellParams[5][4];
    migrationSpeed = cellParams[6][4];
    baseKillProb = cellParams[7][4];
    killProb = baseKillProb;
    infScale = cellParams[8][4];
    influenceRadius = cellParams[9][4];
    migrationBias = cellParams[10][4];
    divProb_base = 0;
    divProb_base = 0;
    deathScale = cellParams[12][4];
    migScale = cellParams[13][4];

    rmax = 1.5*radius*2;

}

void Cell::nk_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    // Inhibition via direct contact, modeling binding of pd-1 and pd-l1.

    if(state != 8){return;}

    double distance = calcDistance(otherX);
    if(distance <= radius+otherRadius){
        // the line below assumes a one to one binding of pd1 and pdl1. You can only bind the minimum num of molecules of either pd1 or pdl1.
        double percent_binding = std::min(otherpdl1, pd1_available) / pd1_expression_level; // this turns it into a fraction.
        double rnd = master_rng.uniform(0,1,temporary_rng);
        if(rnd < percent_binding ){
            next_killProb = next_killProb * percent_binding * inhibitory_effect_of_binding_PD1_PDL1;
            next_migrationSpeed = next_migrationSpeed * percent_binding * inhibitory_effect_of_binding_PD1_PDL1;
            next_death_prob = next_death_prob * (1 + percent_binding) * inhibitory_effect_of_binding_PD1_PDL1;
        }
    }
}



void Cell::nk_update_properties_indirect(size_t step_count){
    // Indirect interactions, models the secretion of stimulatory or inhibitory soluble factors
    // Assumes Cancer cells don't secrete the soluble factors, only other imme cells.
    if(state != 8){return;}

    // posInfluence is M1 + Th
    // negInfluence is M2 + Treg + MDSC
    double posInfluence = 1 - (1 - influences[1])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[5])*(1 - influences[10]);
    double scale = posInfluence - negInfluence;

    next_killProb = next_killProb*pow(infScale, scale);
    next_migrationSpeed = next_migrationSpeed*pow(migScale, scale);
}

void Cell::nk_pd1_expression_level(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug){
    if (type == 4) {
        double temp;
        // Assumes that the expression of PD1 is induced in a contact independent manner, by M2 macs, cancer, Tregs and MDSCs.
        double influence = 1 - (1-influences[2]) * (1-influences[3])*(1-influences[5])*(1-influences[10]);
        if (influence >= threshold_for_pd1_induction ) {
            temp = pd1_expression_level + dt * (influence - threshold_for_pd1_induction)/threshold_for_pd1_induction;
            pd1_expression_level = (temp < max_pd1_level) ? temp : max_pd1_level;
        }
        else {
            temp = pd1_expression_level - pd1_decay_rate * dt;
            pd1_expression_level = (temp > 0) ? temp : 0; // pd1 expression must be non-negative
        }

        double fraction_pd1_bound_by_drug = sensitivity_to_antiPD1(anti_pd1_concentration,binding_rate_pd1_drug);

        pd1_drug_bound = pd1_expression_level * fraction_pd1_bound_by_drug;
        pd1_available = pd1_expression_level * (1-fraction_pd1_bound_by_drug);
    }
}

