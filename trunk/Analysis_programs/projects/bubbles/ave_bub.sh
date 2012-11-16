#!/bin/bash

start_conf=100
end_conf=339
nbub=6
nbubh=$(($nbub/2))

name_eo=(evn odd)

for iconf in $(seq $start_conf $end_conf)
do
    for eo in 0 1
    do
	name_out=0$iconf/bubble_${name_eo[$eo]}
	chris_name_out=0$iconf/chris_bubble_${name_eo[$eo]}
	
	rm -f $chris_name_out
	
        lista_nomi_bubbles=$(for i in $(seq -w $eo 2 $(($nbub-1))); do echo 0$iconf/bubble_000$i; done)
	paste $lista_nomi_bubbles|awk '{if($1=="#"){printf(" ");for(i=0;i<NF/'$nbubh';i++) printf("%s ",$(i+1));printf("\n");}else{a=b=0;for(i=0;i<NF/2;i++){a+=$(i*2+1);b+=$(i*2+2)};if(a!=0)printf("%+016.16lg\t%+016.16lg\n",a/'$nbubh',b/'$nbubh');else printf("\n");}}' |awk '($1=="#"){if(NF>2){hea=$0;printf("%s\n\n",hea) >> "'$chris_name_out'"}else{if(substr($2,1,5)=="CHRIS")ch=51}}{if(ch>0){ch--;print $0 >> "'$chris_name_out'"}else{print $0}}' > $name_out

    done
    
    paste 0$iconf/bubble_evn 0$iconf/bubble_odd|awk -f /Users/francesco/Prace/sunpho/Analysis_programs/projects/bubbles/conv.awk > 0$iconf/bubble_prd
    paste 0$iconf/chris_bubble_evn 0$iconf/chris_bubble_odd|awk -f /Users/francesco/Prace/sunpho/Analysis_programs/projects/bubbles/conv.awk > 0$iconf/chris_bubble_prd
done
