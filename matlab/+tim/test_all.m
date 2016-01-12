A = randn(3, 1000);

tim.differential_entropy_kl(A);
tim.differential_entropy_kl(A, 'k', 2);

tim.differential_entropy_kl_t(A, 100);
tim.differential_entropy_kl_t(A, 100, 'k', 2);
tim.differential_entropy_kl_t(A, 100, 'k', 2, 'filter', 1);
tim.differential_entropy_kl_t(A, 100, 'k', 2, 'filter', 1);
tim.differential_entropy_kl_t(A, 100, 'k', 2, 'filter', ones(1, 11));
tim.differential_entropy_kl_t(A, 100, 'k', 2, 'filter', ones(1, 1001));
tim.differential_entropy_kl_t(A, 10000, 'k', 2, 'filter', 1);

tim.renyi_entropy_lps(A);
tim.renyi_entropy_lps(A, 'q', 2);
tim.renyi_entropy_lps(A, 'q', 2, 'kSuggestion', 2);

tim.renyi_entropy_lps_t(A, 100);
tim.renyi_entropy_lps_t(A, 100, 'q', 2);
tim.renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2);
tim.renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tim.renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tim.renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 11));
tim.renyi_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 1001));
tim.renyi_entropy_lps_t(A, 10000, 'q', 2, 'kSuggestion', 2, 'filter', 1);

tim.tsallis_entropy_lps(A);
tim.tsallis_entropy_lps(A, 'q', 2);
tim.tsallis_entropy_lps(A, 'q', 2, 'kSuggestion', 2);
tim.tsallis_entropy_lps(A, 'q', 2, 'kSuggestion', 2);

tim.tsallis_entropy_lps_t(A, 100);
tim.tsallis_entropy_lps_t(A, 100, 'q', 2);
tim.tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2);
tim.tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tim.tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', 1);
tim.tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 11));
tim.tsallis_entropy_lps_t(A, 100, 'q', 2, 'kSuggestion', 2, 'filter', ones(1, 1001));
tim.tsallis_entropy_lps_t(A, 10000, 'q', 2, 'kSuggestion', 2, 'filter', 1);

tim.differential_entropy_nk(A);

tim.differential_entropy_sp(A);

tim.divergence_wkv(A, A);

tim.mutual_information(A, A);
tim.mutual_information(A, A, 'xLag', 0, 'yLag', 0);
tim.mutual_information(A, A, 'xLag', 0, 'yLag', 0, 'k', 2);

tim.mutual_information_t(A, A, 100);
tim.mutual_information_t(A, A, 100, 'xLag', 0, 'yLag', 0);
tim.mutual_information_t(A, A, 100, 'xLag', 0, 'yLag', 0, 'k', 2);
tim.mutual_information_t(A, A, 100, 'filter', 1);
tim.mutual_information_t(A, A, 100, 'filter', ones(1, 11));
tim.mutual_information_t(A, A, 100, 'filter', ones(1, 1001));
tim.mutual_information_t(A, A, 10000, 'yLag', 2, 'filter', 1);

tim.mutual_information_p(A, A, A);
tim.mutual_information_p(A, A, A, 'xLag', 0, 'yLag', 0, 'zLag', 0);
tim.mutual_information_p(A, A, A, 'k', 2);

tim.mutual_information_pt(A, A, A, 100);
tim.mutual_information_pt(A, A, A, 100, ...
    'xLag', 0, 'yLag', 0, 'zLag', 0);
tim.mutual_information_pt(A, A, A, 100, 'k', 2);
tim.mutual_information_pt(A, A, A, 100, 'k', 2, 'filter', 1);

tim.transfer_entropy(A, A, A);
tim.transfer_entropy(A, A, A, 'xLag', 0, 'yLag', 0, 'wLag', 0);
tim.transfer_entropy(A, A, A, 'xLag', 0, 'yLag', 0, 'wLag', 0, 'k', 2);

tim.transfer_entropy_t(A, A, A, 100);
tim.transfer_entropy_t(A, A, A, 100, 'xLag', 0, 'yLag', 0, 'wLag', 0);
tim.transfer_entropy_t(A, A, A, 100, 'xLag', 0, 'yLag', 0, 'wLag', 0, 'k', 2);
tim.transfer_entropy_t(A, A, A, 100, 'filter', 1);
tim.transfer_entropy_t(A, A, A, 100, 'filter', ones(1, 11));
tim.transfer_entropy_t(A, A, A, 100, 'filter', ones(1, 1001));
tim.transfer_entropy_t(A, A, A, 10000, 'yLag', 2, 'filter', 1);

tim.transfer_entropy_p(A, A, A, A);
tim.transfer_entropy_p(A, A, A, A, 'xLag', 0, 'yLag', 0, 'zLag', 0);
tim.transfer_entropy_p(A, A, A, A, 'k', 2);

tim.transfer_entropy_pt(A, A, A, A, 100);
tim.transfer_entropy_pt(A, A, A, A, 100, ...
    'xLag', 0, 'yLag', 0, 'zLag', 0, 'wLag', 0);
tim.transfer_entropy_pt(A, A, A, A, 100, 'k', 2);
tim.transfer_entropy_pt(A, A, A, A, 100, 'k', 2, 'filter', 1);
