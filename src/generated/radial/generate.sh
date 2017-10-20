sed -i -e 's/x\*\*2/x2/g' radial_gen.part2
sed -i -e 's/x\*\*3/(x2*x)/g' radial_gen.part2
sed -i -e 's/x\*\*4/(x2*x2)/g' radial_gen.part2
sed -i -e  's/x\*\*5/(x2*x2*x)/g' radial_gen.part2
sed -i -e  's/x\*\*6/(x2*x2*x2)/g' radial_gen.part2
sed -i -e  's/x\*\*7/(x2*x2*x2*x)/g' radial_gen.part2
sed -i -e  's/x\*\*8/(x2*x2*x2*x2)/g' radial_gen.part2
sed -i -e  's/x\*\*9/(x2*x2*x2*x2*x)/g' radial_gen.part2
sed -i -e  's/x\*\*10/(x2*x2*x2*x2*x2)/g' radial_gen.part2

sed -i -e  's/y\*\*2/y2/g' radial_gen.part2
sed -i -e  's/y\*\*3/(y2*y)/g' radial_gen.part2
sed -i -e  's/y\*\*4/(y2*y2)/g' radial_gen.part2
sed -i -e  's/y\*\*5/(y2*y2*y)/g' radial_gen.part2
sed -i -e  's/y\*\*6/(y2*y2*y2)/g' radial_gen.part2
sed -i -e  's/y\*\*7/(y2*y2*y2*y)/g' radial_gen.part2
sed -i -e  's/y\*\*8/(y2*y2*y2*y2)/g' radial_gen.part2
sed -i -e  's/y\*\*9/(y2*y2*y2*y2*y)/g' radial_gen.part2
sed -i -e  's/y\*\*10/(y2*y2*y2*y2*y2)/g' radial_gen.part2

sed -i -e  's/p\*\*2/p2/g' radial_gen.part2
sed -i -e  's/p\*\*3/(p2*p)/g' radial_gen.part2
sed -i -e  's/p\*\*4/(p2*p2)/g' radial_gen.part2
sed -i -e  's/p\*\*5/(p2*p2*p)/g' radial_gen.part2
sed -i -e  's/p\*\*6/(p2*p2*p2)/g' radial_gen.part2
sed -i -e  's/p\*\*7/(p2*p2*p2*p)/g' radial_gen.part2
sed -i -e  's/p\*\*8/(p2*p2*p2*p2)/g' radial_gen.part2
sed -i -e  's/p\*\*9/(p2*p2*p2*p2*p)/g' radial_gen.part2
sed -i -e  's/p\*\*10/(p2*p2*p2*p2*p2)/g' radial_gen.part2

rm *-e 
cat radial_gen.part* > radial_gen.cpp
cp radial_gen.cpp ../../lib