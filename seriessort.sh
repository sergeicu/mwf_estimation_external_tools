#!/bin/bash 

dcmtk="/opt/el7/pkgs/dcmtk/3.6.1-20161102/bin/"
outdir="$1"

export dcmtk outdir

get_tag() { $dcmtk/dcmdump "$1" | grep "$2" | head -1 | awk '$0=$2' FS=[ RS=] | tr " " "_"; }

sortd() {
series=`get_tag "$1" 'SeriesDescription'`
seriesnum=`get_tag "$1" 'SeriesNumber'`
dpath="$outdir"/"$seriesnum"_"$series"
if [[ ! "$seriesnum" == "" ]] && [[ ! "$series" == "" ]] ; then
 if [[ ! -d "$dpath" ]]; then
  mkdir -vp "$dpath"
  mv -v --backup=t "$1" "$dpath"
 else
  mv -v --backup=t "$1" "$dpath"
 fi
else
 echo missing dicom 'info', doing nothing 
fi;
}

export -f get_tag sortd

find $outdir -maxdepth 1 -type f -print | /home/ch163210/bin/parallel -j `nproc` -k sortd
