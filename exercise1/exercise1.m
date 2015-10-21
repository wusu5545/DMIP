clear all;
clc;
A=[11 10 14;12 11 -13;14 13 -66];
det(A);
[U S V]=svd(A);
[rows cols]=size(S);
n=min(rows,cols);
Sinv=zeros(size(S));
eps=0.001;

for d = 1:n
    if (S(d,d)>eps)
        Sinv(d,d)=1/S(d,d);
    end
end

Ainv=V*Sinv*U';
inv(A);
cond(A);%returns the 2-norm condition number (the ratio of the
        %largest singular value of X to the smallest).  Large condition
        %numbers indicate a nearly singular matrix.
rank(A);
null(A);
%eigshow

%exercise 1.1
b=[1.001;0.999;1.001];
br=round(b);
xr=Ainv*br;
x=Ainv*b;
xn=(xr-x)./x;

%exercise 1.2
%Optimization 1
nd=[(S(1,1)+S(2,2))/2;(S(1,1)+S(2,2))/2;0];
A0=U*diag(nd)*V';
norm(A,'fro');
norm(A0,'fro');

%Optimization 4
AA=[3 1;2 1;1 1;0 1;-1 1;-1 1;-2 1];
b=[2 1 2 0 1 -1 -1]';
l=pinv(AA)*b;
%[m,t]
figure(1);
plot(AA(:,1),b(:),'.');
hold on;
x=-3:0.1:4;
y=l(1)*x+l(2);
plot(x,y,'.r');
hold off;

%Optimization 2
b=[1 -1 1 -1;1 2 -3 -4];
no=size(b,2);
MM=zeros(no,4);
for i=1:no
    MM(i,1)=b(1,i)^2;
    MM(i,2)=b(1,i)*b(2,i);
    MM(i,3)=MM(i,2);
    MM(i,4)=b(2,i)^2;
end

a=null(MM);
A=reshape(a,2,2)';
norm(A,'fro')
a'*MM*a;

%Optimization 3
xray=imread('yu_fill.jpg');
xray=double(xray);
figure(2);
subplot(1,4,1);
imagesc(xray);
axis image;
title('Original image');
xlabel('x');
ylabel('y');

[U S V]=svd(xray);
I1=U(:,1)*S(1,1)*V(:,1)';
subplot(1,4,2);
colormap gray;
imagesc(I1);
axis image;
title('Rank 1 Approximation');
xlabel('x');
ylabel('y');

subplot(1,4,3);
colormap gray;
imagesc(xray./I1);
axis image;
title('Corrected:xray./I1');
xlabel('x');
ylabel('y');

k=100;
Id=zeros(size(xray));
figure(3);
for d=1:k
    Id=Id+U(:,d)*S(d,d)*V(:,d)';
    colormap gray;
    imagesc(Id);
    axis image;
    title(['Rank ',int2str(d),'Approximation']);
    xlabel('x');
    ylabel('y');
    F(d)=getframe(gcf);
end

movie2avi(F,'rank.avi','fps',5);

figure(2);
subplot(1,4,4)
colormap gray;
imagesc(Id);
axis image;
title(['Rank ','100','Approximation']);
xlabel('x');
ylabel('y');

%Fourier Transform
p=phantom(128);
figure(4);
subplot(2,2,1);
imagesc(p);
axis image;
colormap gray;
title('Phantom Image');

I=fft2(p);
subplot(2,2,2);
imagesc(abs(I));
axis image;

subplot(2,2,3);
imagesc(abs(fftshift(I)));
axis image;
subplot(2,2,4);
imagesc(log(1+abs(fftshift(I))));
axis image;