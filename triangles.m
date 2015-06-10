face1 = [81,80,72];
a = sqrt(face1(1) + face1(2));
b = sqrt(face1(2) + face1(3));
c = sqrt(face1(3) + face1(1));
p = (a + b + c)/2;
area = sqrt(p*(p-a)*(p-b)*(p-c))


