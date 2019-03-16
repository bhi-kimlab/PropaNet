# ../exp/exp.cold.original : expression data
# GroupInfo.cold : group information
# CompareGroup.cold : what group to test (control vs. test)

#Rscript DEGlimma.R ../exp/exp.cold.original GroupInfo.cold CompareGroup.cold signed_zstats ../exp/exp.cold.DEG.zvalue
Rscript DEGlimma.R ../exp/exp.heat.original GroupInfo.heat CompareGroup.heat signed_zstats ../exp/exp.heat.DEG.zvalue

#Rscript DEGlimma.R ../exp/exp.cold.D2 GroupInfo.cold.D2 CompareGroup.cold.D2 signed_zstats ../exp/exp.cold.D2.DEG.zvalue
#Rscript DEGlimma.R ../exp/exp.cold.D3 GroupInfo.cold.D3 CompareGroup.cold signed_zstats ../exp/exp.cold.D3.DEG.zvalue

#Rscript DEGlimma.R ../exp/exp.cold.original GroupInfo.cold CompareGroup.cold signed_binary ../exp/exp.cold.DEG.binary
#Rscript DEGlimma.R ../exp/exp.heat.original GroupInfo.heat CompareGroup.heat signed_binary ../exp/exp.heat.DEG.binary

Rscript DEGlimma.R ../exp/exp.heat.D2 GroupInfo.heat.D2 CompareGroup.heat signed_zstats ../exp/exp.heat.D2.DEG.binary
#Rscript DEGlimma.R ../exp/exp.heat.D3 GroupInfo.heat.D3 CompareGroup.heat signed_binary ../exp/exp.heat.D3.DEG.binary
