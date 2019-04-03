mkdir 10
cd 10
NG=10;
potential=tip4p;
Nsobol=200;
freq=10;
# make directories whole numbers: 1, 20, going to run 20 iterations
# first value is start, and increase by increment
n=1;
max=5;
start=0.0;
inc=0.5;
while [ "$n" -le "$max" ]; do
  mkdir "$n"
  cp ../cage.xyz "$n"
  cp ../submit.sh "$n"
  cd "$n"
#single arrow makes a file, double arrow appends to an existing file
  echo "$potential" > input
#cant use floating point with bash, use bc instead
  echo "$start + $n*$inc" | bc  >> input
  echo "$NG" >> input
  echo "$Nsobol" >> input
  echo "$freq" >> input
  echo cage.xyz >> input
  echo 1 >> input
  cd ..
  n=`expr "$n" + 1`;
done
