clc;
Rx = rotx(0.1);
Ry = roty(0.1);
Rz = rotz(25);

rotm = Rz*Ry*Rx;
eulXYZ = rotm2eul(rotm,'XYZ');
eulZYX = rotm2eul(rotm,'ZYX');

eulXYZindeg = eulXYZ * 180 / pi
eulZYXindeg = eulZYX * 180 / pi