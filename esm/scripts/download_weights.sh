# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

#! /bin/sh
# Usage: bash <repo>/esm/scripts/download_weights.sh /path/to/weights/
# re run the command to continue downloading incomplete files.
echo $0
run_dir=$(readlink -f $(dirname $0)) ;

if [[ "$(which aria2c)" == "" ]];then
  echo aria2c is required for downloading!
  exit 1
fi

for _sh_exec in /bin/sh bash zsh;do
  # only bash and zsh are tested.
  _sh_exec=$(which ${_sh_exec})
  if [[ "$(${_sh_exec} --version)" =~ "bash" || "$(${_sh_exec} --version)" =~ "zsh" ]];then
    echo $(which ${_sh_exec});
    break
  fi
done

model_pth=$1
if [[ "$model_pth" == "" ]];then
  model_pth=$PWD
fi

mkdir -p $model_pth/checkpoints

pushd $model_pth/checkpoints
cat $run_dir/../README.md |grep -e '^|' |grep -e 'fair-esm/models' |tr -d '|' | \
  awk 'BEGIN{print "set -e"};
  {
      # read urls starts with https
      split($0,arr,"https");
      url="https"arr[2];

      # remove blank spaces after url string
      split(url,url_arr," ")
      url=url_arr[1]

      # guessing regression pt url
      url_regression=url;
      sub("models","regression",url_regression);
      sub(".pt","-contact-regression.pt",url_regression);

      # downloading weight
      url_basename_idx=split(url,url_arr,"/")
      url_basename=url_arr[url_basename_idx]
      print "if [[ ! -f "url_basename" || -f "url_basename".aria2 ]];then echo Download not complete: "url_basename";aria2c -x 10 "url";else echo Download complete: "url_basename";fi"

      # downloading regression, if not existing and through an error, we just ignore it.
      url_regression_basename_idx=split(url_regression,url_arr,"/")
      url_regression_basename=url_arr[url_regression_basename_idx]
      print "if [[ ! -f "url_regression_basename" || -f "url_regression_basename".aria2 ]];then echo Download not complete: "url_regression_basename";aria2c -x 10 "url_regression" 2>/dev/null || echo Never mind. "url_regression_basename" may not exist. ;else echo Download complete: "url_regression_basename";fi"

      }' |$_sh_exec


echo "Your model directory is located at \`$(readlink -f ${model_pth})\`."

popd