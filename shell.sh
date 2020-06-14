
../plink_linux_x86_64_20181202/plink --file zqr_new_plink --out zqr_new_plink --allow-extra-chr

../plink_linux_x86_64_20181202/plink --bfile "zqr_new_plink" --r2 --out "zqr_new_plink" --allow-extra-chr



#../mapthin-v1.11-linux-x86_64/mapthin -b 20 zqr_new_plink.bim zqr_new_plink-thin.bim


../plink_linux_x86_64_20181202/plink --bfile "zqr_new_plink" --r2 --ld-window-r2 0.01 --ld-window 999999 --ld-window-kb 8000 --out "zqr_new_plink" --allow-extra-chr --not-chr NW_011499845.1 # r2 0
#../plink_linux_x86_64_20181202/plink --bfile "zqr_new_plink" --r2 --ld-window-r2 0.2 --ld-window 999999 --ld-window-kb 8000 --out "zqr_new_plink" --allow-extra-chr --not-chr NW_011499845.1 # r2 0.2


cat zqr_new_plink.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > zqr_new_plink.ld.summary
