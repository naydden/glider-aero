function tests = TestVelocity
tests = functiontests(localfunctions);
end
% below
function testVelSegment(testCase)
    x_b=0; y_b=0; z_b=0; x_c=0; y_c=1; z_c=0; x=0.5; y=0.5; z=0;
    actSolution = vel_segment (x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
    expSolution = [0;0;-0.225079];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
%above
function testVelSegment2(testCase)
    x_b=0; y_b=0; z_b=0; x_c=0; y_c=1; z_c=0; x=-0.5; y=0.5; z=0;
    actSolution = vel_segment (x_b,y_b,z_b,x_c,y_c,z_c,x,y,z);
    expSolution = [0;0;0.225079];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
% to the right for AB
function testVelImaginaryAB(testCase)
    x_e=0; y_e=0; z_e=0; x=0.5; y=0.5; z=0;
    actSolution = vel_seminfline(x_e,y_e,z_e,x,y,z,'AB');
    expSolution = [0;0;-0.27169];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
%to the left for CD
function testVelImaginaryCD(testCase)
    x_e=0; y_e=1; z_e=0; x=0.5; y=0.5; z=0;
    actSolution = vel_seminfline(x_e,y_e,z_e,x,y,z,'CD');
    expSolution = [0;0;-0.27169];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
% to the left for AB
function testVelImaginaryAB2(testCase)
    x_e=0; y_e=0; z_e=0; x=0.5; y=-0.5; z=0;
    actSolution = vel_seminfline(x_e,y_e,z_e,x,y,z,'AB');
    expSolution = [0;0;0.27169];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
% far to the left of CD
function testVelImaginaryCD2(testCase)
    x_e=0; y_e=1; z_e=0; x=0.5; y=-0.5; z=0;
    actSolution = vel_seminfline(x_e,y_e,z_e,x,y,z,'CD');
    expSolution = [0;0;-0.069828];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end