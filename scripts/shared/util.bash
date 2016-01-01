

test_dir() {
  dir=$1
  label=$2
  if [ ! -d $dir ]; then
    echo "$0 :: error: $label dir does not exist"
    exit
  elif (for f in $(ls $dir); do exit 1; done) then
    echo "$0 :: error: $label dir is empty"
    rm -rf $dir; exit
  elif ! (cd $dir; for f in $(ls); do if ! test -s $f; then exit 1; fi;  done); then
    echo "$0 :: error: $label dir contains empty files"
    rm -rf $dir; exit
  fi
}
