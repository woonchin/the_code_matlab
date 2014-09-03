function [A,B]=a7(N)
M=N/2;
A=zeros(M);
B=zeros(M);
Q=zeros(N,M+1);
A(:,1)=1:2:(N-1);
B(:,1)=2:2:N;
Q(A(:,1),1)=B(:,1);
Q(B(:,1),1)=A(:,1);
S=perms(1:N);
for i=1:M
    S(:,[2*i-1,2*i])=[min(S(:,[2*i-1,2*i]),[],2),max(S(:,[2*i-1,2*i]),[],2)];
end
S=unique(S,'rows');
for s=1:(M-1)  
    for i=1:M
        S(S(:,2*i-1)==A(i,s),:)=[];
        S(S(:,2*i)==A(i,s),:)=[];
        S(S(:,2*i-1)==B(i,s),:)=[];
        S(S(:,2*i)==A(i,s),:)=[];
    end
p=zeros(size(S),1);
for k=1:M
    for i=1:size(S)
         P(i)=0
            for j=1:M
                if S(i,2*j-1)==A(k,s)&& S(i,2*j)==B(k,s)
                    P(i)=1;
                    break;
                end
            end
        end
    end
end
S(P==1,:)=[];
x=S(S(:,2*s+1)==1,:);
r=randi(size(x,1),1);
A(:,s+1)=x(r,1:2:(N-1));
B(:,s+1)=x(r,2:2:N);
A
B
size(x)
          
    
end

