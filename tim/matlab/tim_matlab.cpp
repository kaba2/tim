// Documentation: tim_matlab_cpp.txt

#define FORCE_LINKING(name) \
	void force_linking_##name(); \
	force_linking_##name()

namespace
{

	void forceLinking()
	{
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
	}

}
