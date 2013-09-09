// DocumentationOf: tim_matlab.h

/*
This file is used as the sole file to build against the
TIM libraries. The idea is to minimize the amount of the
build work done on the Matlab side when building the TIM
Matlab interface.

It is necessary for this file to refer to things in
other translation units. This causes the translation
units to self-register themselves to PastelMatlab.
*/

#define FORCE_LINKING(name) \
	void force_linking_##name(); \
	void call_force_linking_##name() \
	{ \
		force_linking_##name(); \
	}

FORCE_LINKING(differential_entropy_kl);
FORCE_LINKING(differential_entropy_kl_t);
FORCE_LINKING(differential_entropy_nk);
FORCE_LINKING(differential_entropy_sp);
FORCE_LINKING(divergence_wkv);
FORCE_LINKING(entropy_combination);
FORCE_LINKING(entropy_combination_t);
FORCE_LINKING(mutual_information_naive);
FORCE_LINKING(renyi_entropy_lps);
FORCE_LINKING(renyi_entropy_lps_t);
FORCE_LINKING(tsallis_entropy_lps);
FORCE_LINKING(tsallis_entropy_lps_t);
