

cat foo | awk '{a=$1;b=$2;c=$3; if(a>b) {x=a; a=b; b=x;} if(b>c) {x=b; b=c; c=x;} if(a>b) {x=a; a=b; b=x;} print a,b,c;}' | sort -k1n | uniq | wc
