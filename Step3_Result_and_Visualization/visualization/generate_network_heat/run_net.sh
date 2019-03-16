
out=$1 #cold/heat
tpCnt=$2 #The number of timepoint

#for ((i=1;i<=$tpCnt;i++));do 
#cut -f1 "../exp.$out.D2.DEG.zvalue.$i.trim"|awk 'NR>1 {print}'>"$out.D2.t$i.TF"; done

#for ((i=1;i<=$tpCnt;i++));do 
#awk 'FNR==1 && NR!=1 {print}{print}' ../$out.D2.t$i.*.TG|sort|uniq > "$out.D2.t$i.TG"; done

#for ((i=1;i<=$tpCnt;i++));do 
#awk 'FNR==1 && NR!=1 {print}{print}' "$out.D2.t$i.*"|sort|uniq > "concat.$out.D2.t$i"; done

for((i=1;i<=$tpCnt;i++));
do
python finalNetwork.py "../../network/templateNetwork.DEG.$out.D2.$i" -bin "exp.$out.DEG.binary" -TF "$out.D2.t$i.TF" "concat.$out.D2.t$i" "$out.D2" -n "$i"
done

python finalNetwork.py "../../network/templateNetwork.$out.D2.corr" "concat.$out.D2" "$out.D2"


