#!/bin/bash

sed 's/\[/(/g;s/\]/)/g;s/Sin/sin/g;s/Cos/cos/g;s/1,/1.,/g;s/2,/2.,/g;s/3,/3.,/g;s/4,/4.,/g;s/8,/8.,/g;s/16,/16.,/g;s/1)/1.)/g;s/2)/2.)/g;s/4)/4.)/g;s/8)/8.)/g;;s/16)/16.)/g;s/exx/e_xx/g;s/eyy/e_yy/g;s/eyx/e_yx/g;s/e_yx\^2/e_yx\*e_yx/g;s/ //g' ./temp.txt > ./out.txt

sed -e :x -e 'N;s/\n//;tx' ./out.txt > ./ot.txt


