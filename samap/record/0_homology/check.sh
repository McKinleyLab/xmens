hu=$(grep -ch "^# CPU time" /n/eddy_lab/users/nhilgert/data/proteomes/hm/hu_vs_mu/*.log 2>/dev/null | paste -sd+ | bc)
mu=$(grep -ch "^# CPU time" /n/eddy_lab/users/nhilgert/data/proteomes/hm/mu_vs_hu/*.log 2>/dev/null | paste -sd+ | bc)
running=$(squeue -u $USER -h -t R | wc -l)
printf "hu_vs_mu: %d / 123845 (%d%%)\n" "$hu" "$((hu*100/123845))"
printf "mu_vs_hu: %d / 66558  (%d%%)\n"  "$mu" "$((mu*100/66558))"
printf "Tasks running: %d / 64\n" "$running"
