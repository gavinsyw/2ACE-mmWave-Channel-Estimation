function [Evaluation, rss_eval] = Evaluate_rss(H, cb_test, rss_test)

rss_eval = abs(cb_test * H);
rss_error = abs(rss_eval - rss_test);

Evaluation = mean(rss_error ./ rss_test);

end