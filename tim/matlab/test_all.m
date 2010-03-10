A = randn(3, 1000);

differential_entropy_kl(A);
differential_entropy_kl(A, 2);
differential_entropy_kl(A, 2, 8);

differential_entropy_kl_t(A, 100);
differential_entropy_kl_t(A, 100, 2);
differential_entropy_kl_t(A, 100, 2, [1]);
differential_entropy_kl_t(A, 100, 2, [1], 8);
differential_entropy_kl_t(A, 100, 2, ones(1, 11), 8);
differential_entropy_kl_t(A, 100, 2, ones(1, 1001), 8);
differential_entropy_kl_t(A, 10000, 2, [1], 8);

renyi_entropy_lps(A);
renyi_entropy_lps(A, 2);
renyi_entropy_lps(A, 2, 2);
renyi_entropy_lps(A, 2, 2, 8);

renyi_entropy_lps_t(A, 100);
renyi_entropy_lps_t(A, 100, 2);
renyi_entropy_lps_t(A, 100, 2, 2);
renyi_entropy_lps_t(A, 100, 2, 2, [1]);
renyi_entropy_lps_t(A, 100, 2, 2, [1], 8);
renyi_entropy_lps_t(A, 100, 2, 2, ones(1, 11), 8);
renyi_entropy_lps_t(A, 100, 2, 2, ones(1, 1001), 8);
renyi_entropy_lps_t(A, 10000, 2, 2, [1], 8);

tsallis_entropy_lps(A);
tsallis_entropy_lps(A, 2);
tsallis_entropy_lps(A, 2, 2);
tsallis_entropy_lps(A, 2, 2, 8);

tsallis_entropy_lps_t(A, 100);
tsallis_entropy_lps_t(A, 100, 2);
tsallis_entropy_lps_t(A, 100, 2, 2);
tsallis_entropy_lps_t(A, 100, 2, 2, [1]);
tsallis_entropy_lps_t(A, 100, 2, 2, [1], 8);
tsallis_entropy_lps_t(A, 100, 2, 2, ones(1, 11), 8);
tsallis_entropy_lps_t(A, 100, 2, 2, ones(1, 1001), 8);
tsallis_entropy_lps_t(A, 10000, 2, 2, [1], 8);

differential_entropy_nk(A);
differential_entropy_nk(A, 8);

differential_entropy_sp(A);

divergence_wkv(A, A);
divergence_wkv(A, A, 8);

mutual_information(A, A);
mutual_information(A, A, 0, 0);
mutual_information(A, A, 0, 0, 2);
mutual_information(A, A, 0, 0, 2, 8);

mutual_information_t(A, A, 100);
mutual_information_t(A, A, 100, 0, 0);
mutual_information_t(A, A, 100, 0, 0, 2);
mutual_information_t(A, A, 100, 0, 0, 2, [1]);
mutual_information_t(A, A, 100, 0, 0, 2, ones(1, 11), 8);
mutual_information_t(A, A, 100, 0, 0, 2, ones(1, 1001), 8);
mutual_information_t(A, A, 10000, 0, 0, 2, [1], 8);

mutual_information_p(A, A, A);
mutual_information_p(A, A, A, 0, 0, 0);
mutual_information_p(A, A, A, 0, 0, 0, 2);
mutual_information_p(A, A, A, 0, 0, 0, 2, 8);

mutual_information_pt(A, A, A, 100);
mutual_information_pt(A, A, A, 100, 0, 0, 0);
mutual_information_pt(A, A, A, 100, 0, 0, 0, 2);
mutual_information_pt(A, A, A, 100, 0, 0, 0, 2, [1]);
mutual_information_pt(A, A, A, 100, 0, 0, 0, 2, [1], 8);

