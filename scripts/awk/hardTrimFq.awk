(FNR-1) % 2 == 0 { name=$1; tag=$2; next }
(FNR-2) % 4 == 0 { substr($0,1,s) substr($0,S) }
                 { print name tag
                   print substr($0,1,s) substr($0,S,E) }
