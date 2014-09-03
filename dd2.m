Y=999;
Z=2999;
K=4999;
P=6999;
for t=1:N/2
    a1(t)=vec_Ez(2*t-1,Y);
    a2(t)=vec_Ez(2*t-1,Y+1);
    a3(t)=vec_Ez(2*t-1,Y+2);
    a4(t)=vec_Ez(2*t-1,Y+3);
    a5(t)=vec_Ez(2*t-1,Y+4);

    
    b1(t)=imag(vec_Hy(2*t,Y));
    b2(t)=imag(vec_Hy(2*t,Y+1));
    b3(t)=imag(vec_Hy(2*t,Y+2));
    b4(t)=imag(vec_Hy(2*t,Y+3));
    b5(t)=imag(vec_Hy(2*t,Y+4));
end

subplot(5,1,1);
plot(a1);
title('Ez-999');
subplot(5,1,2);
plot(a2);
subplot(5,1,3);
plot(a3);
subplot(5,1,4);
plot(a4);
subplot(5,1,5);
plot(a5);
figure;

subplot(5,1,1);
plot(b1);
title('Hy-999');
subplot(5,1,2);
plot(b2);
subplot(5,1,3);
plot(b3);
subplot(5,1,4);
plot(b4);
subplot(5,1,5);
plot(b5);
figure;

for t=1:N/2
    a1(t)=vec_Ez(2*t-1,Z);
    a2(t)=vec_Ez(2*t-1,Z+1);
    a3(t)=vec_Ez(2*t-1,Z+2);
    a4(t)=vec_Ez(2*t-1,Z+3);
    a5(t)=vec_Ez(2*t-1,Z+4);

    
    b1(t)=imag(vec_Hy(2*t,Z));
    b2(t)=imag(vec_Hy(2*t,Z+1));
    b3(t)=imag(vec_Hy(2*t,Z+2));
    b4(t)=imag(vec_Hy(2*t,Z+3));
    b5(t)=imag(vec_Hy(2*t,Z+4));
end

subplot(5,1,1);
plot(a1);
title('Ez-2999');
subplot(5,1,2);
plot(a2);
subplot(5,1,3);
plot(a3);
subplot(5,1,4);
plot(a4);
subplot(5,1,5);
plot(a5);
figure;

subplot(5,1,1);
plot(b1);
title('Hy-2999');
subplot(5,1,2);
plot(b2);
subplot(5,1,3);
plot(b3);
subplot(5,1,4);
plot(b4);
subplot(5,1,5);
plot(b5);
figure;

for t=1:N/2
    a1(t)=vec_Ez(2*t,K);
    a2(t)=vec_Ez(2*t,K+1);
    a3(t)=vec_Ez(2*t,K+2);
    a4(t)=vec_Ez(2*t,K+3);
    a5(t)=vec_Ez(2*t,K+4);

    
    b1(t)=imag(vec_Hy(2*t-1,K));
    b2(t)=imag(vec_Hy(2*t-1,K+1));
    b3(t)=imag(vec_Hy(2*t-1,K+2));
    b4(t)=imag(vec_Hy(2*t-1,K+3));
    b5(t)=imag(vec_Hy(2*t-1,K+4));
end

subplot(5,1,1);
plot(a1);
title('Ez-4999');
subplot(5,1,2);
plot(a2);
subplot(5,1,3);
plot(a3);
subplot(5,1,4);
plot(a4);
subplot(5,1,5);
plot(a5);
figure;

subplot(5,1,1);
plot(b1);
title('Hy-4999');
subplot(5,1,2);
plot(b2);
subplot(5,1,3);
plot(b3);
subplot(5,1,4);
plot(b4);
subplot(5,1,5);
plot(b5);
figure;

for t=1:N/2
    a1(t)=vec_Ez(2*t,P);
    a2(t)=vec_Ez(2*t,P+1);
    a3(t)=vec_Ez(2*t,P+2);
    a4(t)=vec_Ez(2*t,P+3);
    a5(t)=vec_Ez(2*t,P+4);

    
    b1(t)=imag(vec_Hy(2*t-1,P));
    b2(t)=imag(vec_Hy(2*t-1,P+1));
    b3(t)=imag(vec_Hy(2*t-1,P+2));
    b4(t)=imag(vec_Hy(2*t-1,P+3));
    b5(t)=imag(vec_Hy(2*t-1,P+4));
end

subplot(5,1,1);
plot(a1);
title('Ez-6999');
subplot(5,1,2);
plot(a2);
subplot(5,1,3);
plot(a3);
subplot(5,1,4);
plot(a4);
subplot(5,1,5);
plot(a5);
figure;

subplot(5,1,1);
plot(b1);
title('Hy-6999');
subplot(5,1,2);
plot(b2);
subplot(5,1,3);
plot(b3);
subplot(5,1,4);
plot(b4);
subplot(5,1,5);
plot(b5);
figure;