#!/bin/bash

# <Replace with the description and/or purpose of this shell script.>

# GLOBAL_VAR1="one"
# GLOBAL_VAR2="two"

# function function_one() {
#   local LOCAL_VAR1="one"
#   # <Replace with function code.>
# }

# Main body of the shell script starts here.
 python network_weight.py -nwk data/templateNetwork -exp data/DEG.AtGenExpress.signed_zstats.heat_shoots -o data/templateNetwork.heat_shoots

python TF_adding_NP_noCtrl.py data/Ath_TF_list.gene data/templateNetwork.heat_shoots data/DEG.AtGenExpress.signed_zstats.heat_shoots data/DEG.AtGenExpress.signed_binary.heat_shoots -cond AtGenExpress.heat_shoots -outD result

# <Replace with the main commands of your shell script.>

# Exit with an explicit exit status.
exit 0

