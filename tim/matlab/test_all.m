A = randn(3, 1000);

differential_entropy_kl(A);
differential_entropy_kl(A, 'k', 2);

differential_entropy_kl_t(A, 100);
differential_entropy_kl_t(A, 100, 'k', 2);
differential_entropy_kl_t(A, 100, 'k', 2, 'filter', 1);
differential_entropy_kl_t(A, 100, 'k', 2, 'filter', 1);
differential_entropy_kl_t(A, 100, 'k', 2, 'filter', ones(1, 11));
differential_entropy_kl_t(A, 100, 'k', 2, 'filter', ones(1, 1001));
differential_entropy_kl_t(A, 10000, 'k', 2, 'filter', 1);

renyi_entropy_lps(A);
renyi_entropy_lps(A, 'q', 2);
renyi_entropy_lps(A, 'q', 2, 'kSuggestion', 2);

renyi_entropy_lps_t(A, 100);
renyi_entropy_lps_t(A, 100, 'q', 2);
renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2);
renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 11));
renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 1001));
renyi_entropy_lps_t(A, 10000, 'q', 2, 'kSuggestion', 2, 'filter', 1);

tsallis_entropy_lps(A);
tsallis_entropy_lps(A, 'q', 2);
tsallis_entropy_lps(A, 'q', 2, 'kSuggestion', 2);
tsallis_entropy_lps(A, 'q', 2, 'kSuggestion', 2);

tsallis_entropy_lps_t(A, 100);
tsallis_entropy_lps_t(A, 100, 'q', 2);
tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2);
tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 11));
tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 1001));
tsallis_entropy_lps_t(A, 10000, 'q', 2, 'kSuggestion', 2, 'filter', 1);

differential_entropy_nk(A);

differential_entropy_sp(A);

divergence_wkv(A, A);

mutual_information(A, A);
mutual_information(A, A, 'xLag', 0, 'yLag', 0);
mutual_information(A, A, 'xLag', 0, 'yLag', 0, 'k', 2);

mutual_information_t(A, A, 100);
mutual_information_t(A, A, 100, 'xLag', 0, 'yLag', 0);
mutual_information_t(A, A, 100, 'xLag', 0, 'yLag', 0, 'k', 2);
mutual_information_t(A, A, 100, 'filter', 1);
mutual_information_t(A, A, 100, 'filter', ones(1, 11));
mutual_information_t(A, A, 100, 'filter', ones(1, 1001));
mutual_information_t(A, A, 10000, 'zLag', 2, 'filter', 1);

mutual_information_p(A, A, A);
mutual_information_p(A, A, A, 'xLag', 0, 'yLag', 0, 'zLag', 0);
mutual_information_p(A, A, A, 'k', 2);

mutual_information_pt(A, A, A, 100);
mutual_information_pt(A, A, A, 100, ...
    'xLag', 0, 'yLag', 0, 'zLag', 0);
mutual_information_pt(A, A, A, 100, 'k', 2);
mutual_information_pt(A, A, A, 100, 'k', 2, 'filter', 1);

transfer_entropy(A, A, A);
transfer_entropy(A, A, A, 'xLag', 0, 'yLag', 0, 'wLag', 0);
transfer_entropy(A, A, A, 'xLag', 0, 'yLag', 0, 'wLag', 0, 'k', 2);

transfer_entropy_t(A, A, A, 100);
transfer_entropy_t(A, A, A, 100, 'xLag', 0, 'yLag', 0, 'wLag', 0);
transfer_entropy_t(A, A, A, 100, 'xLag', 0, 'yLag', 0, 'wLag', 0, 'k', 2);
transfer_entropy_t(A, A, A, 100, 'filter', 1);
transfer_entropy_t(A, A, A, 100, 'filter', ones(1, 11));
transfer_entropy_t(A, A, A, 100, 'filter', ones(1, 1001));
transfer_entropy_t(A, A, A, 10000, 'wLag', 2, 'filter', 1);

transfer_entropy_p(A, A, A, A);
transfer_entropy_p(A, A, A, A, 'xLag', 0, 'yLag', 0, 'zLag', 0);
transfer_entropy_p(A, A, A, A, 'k', 2);

transfer_entropy_pt(A, A, A, A, 100);
transfer_entropy_pt(A, A, A, A, 100, ...
    'xLag', 0, 'yLag', 0, 'zLag', 0, 'wLag', 0);
transfer_entropy_pt(A, A, A, A, 100, 'k', 2);
transfer_entropy_pt(A, A, A, A, 100, 'k', 2, 'filter', 1);
