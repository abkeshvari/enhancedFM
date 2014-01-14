function [T]=FM_ver1(X,m1,m2)
pres=-10;
X=[eye(m1) zeros(m1,m2); zeros(m2,m1) -eye(m2);X]; %Make the matrix
n= size(X,1);% Number of units
x=X(:,1:m1);y=X(:,m1+1:end);
U=[y -x]'; % Matrix U
T=eye(m1+m2);%Intitialization of T
tau=size(T,1);
O=zeros(n,1); % Keeps the order of processed units
D=zeros(n,1); % The vector D is the list of processed and inefficient units
ind=0; % Index of the current unit
Z=T(:,1:m1+m2)*U(:,:);
while sum(D==0)>0 % If there is units which are not processed yet...
    ind=ind+1;
    Unprocessed=find(D==0); % The list of unprocessed units
    if ind>m1+m2
        Z1=Z;
        Z1(Z1<0)=100000000;
        Z1=min(Z1 ,[],1);
        [~,k]=max(Z1);
    else
        k=ind;
    end;
    k=Unprocessed(k); % The unit with maximum efficiency is chosen to be processed
    if ind<=m1+m2 % For the artificial units just follow the order of units
        k=ind;
    else % Remove the inefficient units from the list
        D(Unprocessed(max(Z,[],1)<10^-20))=-1;
    end;
    T(tau+1,size(T,2)+1)=1;
    D(k)=1;
    O(ind)=k;
    tau=size(T,1);
    u=[U(1:m1+m2,k);zeros(ind,1)]; u(end)=1;
    w=T*u;
    d=Normalization(reshape(w,1,tau),tau);
    np=size(w(w>0),1);
    nn=size(w(w<0),1);
    nz=size(w(w==0),1);
    if np==1 && ind>m1+m2
        T(1:tau,m1+m2+ind)=ones(1,tau)./d(1:tau);
        T(tau,:)=[];
        TT=T;
    end;
    if np> 1 || ind<=m1+m2
        % uncomment the following line to see the progress of the algorithm
%         fprintf('iteration: %d, unit: %d, and %d units left \n',ind, k ,(sum(D==0)))

        rows_T=np*nn+nz; %Number of rows of A, #columns is tau
        Z=T(:,1:m1+m2)*U(:,O(1:ind-1 ));
        Z=abs(Z)<10^-10;
        z_p=(Z(w>0,:));z_n= single(Z(w<0,:));
        if ind<=m1+m2
            z=z_p*z_n';
            weak=single(roundn(T(:,1:m1+m2),pres)==0);
            weak=weak(w>0,:)*weak(w<0,:)';
            dimensions=(z+ weak)>=m1+m2-2;
            [row, col]=find(dimensions);
            pairs=[col row];
        else
            rows_z_n=find(sum(z_n(:,sum(z_p(1:end-1,:),1)>0),2)>=m1+m2 );
            z=single(z_p(1:end-1,:)*z_n(rows_z_n,:)');
            dimensions=z>=m1+m2-2;
            [row, col]=find(dimensions);
            if size(row,1)==1
                row=row';
                col=col';
            end;
            pairs1=[rows_z_n(col) row];
            pairs2=[];col=[];row=[];
            rows_z_n=find(sum(z_n(:,sum(z_p(1:end-1,:),1)>0),2)==m1+m2-2 );
            if ~isempty(rows_z_n)
                %to avoid multiplication of very big matrices, we check the
                %size of the matrices, and if it is too large, we get two
                %blockes of matrices and do the multiplications seperately
                if rows_T>100000000
                    block=floor(size(z_p,1)/2);
                    %first block
                    z=single(z_p(1:block,:)*z_n(rows_z_n,:)');
                    dimensions=z>=m1+m2-2;
                    [row, col]=find(dimensions);
                    %second block
                    z=single(z_p(block:end-1,:)*z_n(rows_z_n,:)');
                    dimensions=z>=m1+m2-2;
                    [row1, col1]=find(dimensions);
                    %combining the results of two blocks
                    row=[row;row1+block];
                    col=[col;col1];
                else
                    z=single(z_p(1:end-1,:)*z_n(rows_z_n,:)'); dimensions=((z))>=m1+m2-2;[row, col]=find(dimensions);
                end;
                if size(row  ,1)==1
                    row=row';
                    col=col';
                end;
                pairs2=[rows_z_n(col) row];
            end;
            pairs1=[pairs1;pairs2];
            pairs2=[];col=[];row=[];
            
            rows_z_n=find(sum(z_n(:,sum(z_p(1:end-1,:),1)>0),2)==m1+m2-1 );
            if ~isempty(rows_z_n)
                if rows_T>=100000000
                    block=floor(size(z_p,1)/2);
                    z=single(z_p(1:block,:)*z_n(rows_z_n,:)');
                    dimensions=((z))>=m1+m2-2;
                    [row, col]=find(dimensions);
                    z=single(z_p(block:end-1,:)*z_n(rows_z_n,:)');
                    dimensions=((z))>=m1+m2-2;
                    [row1, col1]=find(dimensions);
                    row=[row;row1+block];
                    col=[col;col1];
                else
                    z=single(z_p(1:end-1,:)*z_n(rows_z_n,:)');
                    dimensions=z>=m1+m2-2;
                    [row, col]=find(dimensions);
                end;
                if size(row  ,1)==1
                    row=row';
                    col=col';
                end;
                pairs2=[rows_z_n(col) row];
            end;
            pairs=[pairs1;pairs2 ;(1:nn)' ones(nn,1)*np];
        end;
        TD=bsxfun(@times,d',T);
        w_p=(1:size(w,1))'.*(w>0);
        w_p(w_p==0)=[];
        w_n=(1:size(w,1))'.*(w<0);
        w_n(w_n==0)=[];
        w_z=(1:size(w,1))'.*(w==0);
        w_z(w_z==0)=[];
        i1=w_n(pairs(:,1));i2=w_p(pairs(:,2));
        T= TD(i1,:)+TD(i2,:);
        T=[T;TD(w_z,:);];
        if ind>m1+m2
            norms = normr(T);
            [~,ia,~]=unique(roundn(norms(:,1:m1+m2),pres),'rows');
            ia=sort(ia);
            T=norms(ia,:);
        end;
 
    end;
    Z=T(:,1:m1+m2)*U(:,D==0);
    tau=size(T,1);
end;
end