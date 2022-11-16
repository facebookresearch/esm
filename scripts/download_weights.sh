#! /bin/sh
# Usage: bash <repo>/esm/scripts/download_weights.sh /path/to/weights/
# re run the command to continue downloading incomplete files.

run_dir=$(readlink -f $(dirname $0)) ;

if [[ "$(which aria2c)" == "" ]];then
  echo aria2c is required for downloading!
  exit 1
fi

model_pth=$1
if [[ "$model_pth" == "" ]];then
  model_pth=$PWD
fi

mkdir -p $model_pth/checkpoints

pushd $model_pth/checkpoints
cat $run_dir/../README.md |grep -e '^|' |grep -e 'fair-esm/models' |tr -d '|' | \
  awk '{
            # read urls started with https
            split($0,arr,"https");
            url="https"arr[2];
            print url;

            # inspect regression pt url
            url_regression=url;
            sub("models","regression",url_regression);
            sub(".pt","-contact-regression.pt",url_regression);
            print url_regression;

            # download weight
            system("if [[ ! -f $(basename "url") || -f $(basename "url").aria2 ]];then echo Download not complete: $(basename "url");aria2c -x 10 "url";else echo Download complete: $(basename "url");fi")

            # download regression
            system("if [[ ! -f $(basename "url_regression") || -f $(basename "url_regression").aria2 ]];then echo Download not complete: $(basename "url_regression");aria2c -x 10 "url_regression";else echo Download complete: $(basename "url_regression");fi")

            }'

popd