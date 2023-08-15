   # add_bt_tag.awk

   BEGIN {
       OFS = "\t"
   }

   {
       if (\$3 != "*" && \$3 != "=")
           \$12 = \$12 "\tBT:Z:" \$3
       else
           \$12 = \$12 "\tBT:Z:NA"
       print
   }