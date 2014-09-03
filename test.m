clear all;
m_x=9; %matrix column size
m_y=9;
n_x=28; %matrix row size
n_y=28;
x_dis_radius=.5;
x_dis_angle=pi/12;
y_dis_radius=.5;
y_dis_angle=pi/12;


X=zeros(m_x,n_x);
Y=zeros(m_y,n_y);

for i=1:m_x
    for j=1:n_y
        X(i,j)=(2+(i-1)*x_dis_radius)*cos((j-1)*x_dis_angle);
        Y(i,j)=(2+(i-1)*y_dis_radius)*sin((j-1)*y_dis_angle);
    end
end

Z=real ( sqrt(4-(sqrt(X.^2+Y.^2)-4).^2));
mesh(X,Y,Z);
hold on;
mesh(X,Y,-Z)
axis equal;