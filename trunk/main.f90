program main

use nr,only:indexx,jacobi,gaussj
      
implicit none

integer,parameter :: n_pop=1000,L=300,n_gen=20000,run_thresh=15,n_gates=20,n_gates_nom=12,n_transient=5,n_inp=4
real,parameter :: p_wire=1.0
!integer,parameter :: n_pop=1000,L=300,n_gen=100000,run_thresh=6,n_gates=50,n_gates_nom=30,n_transient=15,n_inp=8
integer :: goal,gen,gate1,gate2
integer :: a,b,c,run=0,max_qual=0
integer,dimension(n_inp , 2**n_inp) :: inp = 0
integer,dimension(n_gates) :: active = 0
integer,dimension(n_gates,n_gates) :: adj = 0
integer,dimension(n_gates,n_pop) :: x=0,y=0
integer,dimension( 2**n_inp) :: out = 0
integer,dimension(8) :: datevals
integer,dimension( n_pop ) :: order = 0
integer,dimension(n_gates , n_pop):: in1=0 , in2=0
integer,dimension(n_gates , n_pop) :: state=0,oldstate=0
real, dimension(n_pop) :: qual = 0.0,wire_list=0.0,mod_list=0.0
real,dimension(4) :: rand
logical :: again = .true.
integer :: old_sum,epoch=200
real :: b_wire=0.0

character(len=8) ::  date
character(len=10) ::  time
character(len=10) ::  zone
character(len=100) :: filename

call date_and_time(date, time, zone,datevals)
filename = "modularity-"//date//"-"//time(1:6)//".dat"

!CALL RANDOM_SEED
print*,datevals
!CALL RANDOM_SEED 
call system_clock(count=a)
print *,a
CALL RANDOM_SEED(put = datevals*a)
!call random_number(rand)



do a = 1,2**n_inp
   do b = 1,n_inp
      inp(b,a)=transfer(btest(a-1,b),0)
   end do
end do

!out = iand(ior(ieor(inp(1,:),inp(2,:)),ieor(inp(3,:),inp(4,:))),ior(ieor(inp(5,:),inp(6,:)),ieor(inp(7,:),inp(8,:))))
!out = iand(ieor(inp(1,:),inp(2,:)),ieor(inp(3,:),inp(4,:)))
out = ieor(ior(inp(1,:),inp(2,:)),iand(inp(3,:),inp(4,:)))
!out(5) = 0
!print*, out 
!pause
!generate nets
do a = 1,n_pop
again = .true.
do while (again)
   do b = 1,n_inp
      in1(b,a) = b
      in2(b,a) = b
   end do
   
   do b = n_inp+1,n_gates
      call random_number(rand)
          in1(b,a) = floor(rand(1)*n_gates)+1
      in2(b,a) = floor(rand(2)*n_gates)+1 
   end do

   active = 0
   old_sum = 0
   active(n_gates) = 1
   
   do while (sum(active).gt.old_sum)
      old_sum = sum(active)
          do b = 1,n_gates
             if(active(n_gates+1-b).gt.0) then
                    active(in1(n_gates+1-b,a)) = 1
                        active(in2(n_gates+1-b,a)) = 1
                  end if
      end do
   end do
   
   if(sum(active(1:n_inp)) < n_inp) then
      again = .true.
   else
      if (sum(active(n_inp+1:n_gates-1)) /= n_gates_nom) then
             again = .true.
      else
         again = .false.
      end if
   end if 

end do
!end net generation loop
end do

print*,"Go!"

!start main loop
do gen = 1,n_gen


!print*,gen
qual = 0

!evaluate nets
do c = 1,2**n_inp
    state = 0
    oldstate = 0
        x = 0
        y = 0
        do a = 1,n_transient
                state(1:n_inp,:) = spread(inp(:,c),ncopies=n_pop,dim=2)
                oldstate = state
                do b = 1,n_pop
                        x(n_inp+1:n_gates,b) = state(in1(n_inp+1:n_gates,b),b)
                        y(n_inp+1:n_gates,b) = state(in2(n_inp+1:n_gates,b),b)
                end do
                state(n_inp+1:n_gates,:) = 1 - x(n_inp+1:n_gates,:)*y(n_inp+1:n_gates,:)
        end do
        
        !print*,sum(1-(IEOR(state,oldstate)),dim=1)/4

        qual = qual + (1 - (IEOR(state(n_gates,:),spread(out(c),ncopies=n_pop,dim=1))))*sum(1-(IEOR(state,oldstate)),dim=1)/n_gates
                
end do
!print*,qual
!pause
!disqualify nets with wrong number of gates, unconnected inputs

do a = 1,n_pop
   active = 0
   old_sum = 0
   active(n_gates) = 1
   
   do while (sum(active).gt.old_sum)
      old_sum = sum(active)
          do b = 1,n_gates
             if(active(n_gates+1-b).gt.0) then
                 active(in1(n_gates+1-b,a)) = 1
                 active(in2(n_gates+1-b,a)) = 1
             end if
      end do
   end do
   
   if(sum(active(1:n_inp)).lt.n_inp) then
      qual(a) = 0.0
          
   else
      if (sum(active(n_inp+1:n_gates-1)) /= n_gates_nom) then
             qual(a) = 0.0
      end if
   end if 
 max_qual = maxval(qual)
 !end do
 !do a = 1,n_pop

   adj = 0
   do b = 1,n_gates
     adj(b,in1(b,a)) = 1
     adj(b,in2(b,a)) = 1
     adj(in1(b,a),b) = 1
     adj(in2(b,a),b) = 1
     
   end do

   if ((b_wire>0).and.(qual(a)>0) ) then 
           !print *,wire(adj,active)
           wire_list(a) = wire(adj,active)
           
           call random_number(rand)
           if (rand(1)<p_wire) then 
               qual(a) = qual(a) - b_wire*wire_list(a)
           endif
   endif
  
end do
!take measurements at regular intervals


if (mod(gen,epoch)==1) then
    !measure modularity
   
    do a = 1,n_pop
    active = 0
    old_sum = 0
    active(n_gates) = 1
   
    do while (sum(active).gt.old_sum)
      old_sum = sum(active)
          do b = 1,n_gates
             if(active(n_gates+1-b).gt.0) then
                 active(in1(n_gates+1-b,a)) = 1
                 active(in2(n_gates+1-b,a)) = 1
             end if
      end do
    end do
    
    adj = 0
      do b = 1,n_gates
         adj(b,in1(b,a)) = 1
         adj(b,in2(b,a)) = 1
         adj(in1(b,a),b) = 1
         adj(in2(b,a),b) = 1
       end do

    if ((qual(a)>0)) then 
       wire_list(a) = wire(adj,active)
    else
       wire_list(a) = 1
    endif    
    mod_list(a) = modularity(adj)
    end do


    !change environment on certain conditions
    if(maxval(qual)/2.**n_inp >= .99) then
        epoch = 20
        run = run + 1
        !print *,run
        if (run>run_thresh) then
           b_wire = 0.1
        endif
        if (run>(run_thresh*2)) then
           b_wire = 0.0
        endif
        if (run>(run_thresh*3)) then
           exit
        endif
    endif


    print*,gen-1,maxval(qual)/2.**n_inp,real(sum(qual))/(2.**n_inp*n_pop),sum(wire_list)/n_pop, &
           sum(mod_list)/n_pop,maxval(mod_list)
    OPEN(1, FILE=filename,position='append')  
    WRITE (1,*) gen-1,maxval(qual)/2.**n_inp,real(sum(qual))/(2.**n_inp*n_pop),sum(wire_list)/n_pop, & 
                sum(mod_list)/n_pop,maxval(mod_list)
close(1)
endif

!sort nets by quality

call indexx(qual,order)
!print*,maxval(qual)/2.**n_inp,real(sum(qual))/(2.**n_inp*n_pop),sum(wire_list)/n_pop,sum(mod_list)/n_pop

!pause

in1(:,:) = in1(:,order)
in2(:,:) = in2(:,order)
!print*,qual
!qual = qual(order)

!replace worst nets
in1(:,1:L) = in1(:,n_pop-L+1:n_pop)
in2(:,1:L) = in2(:,n_pop-L+1:n_pop)



!mutate nets
do a = 1,n_pop-L
        call random_number(rand)
        gate1 = floor(rand(1)*(n_gates-n_inp))+1+n_inp
        gate2 = floor(rand(2)*(n_gates))+1
        
        if (rand(3).lt.0.5) then
                in1(gate1,a) = gate2
        else
                in2(gate1,a) = gate2
        end if
end do

!end main loop
end do

contains 

function diag(a) result(diagonal)
      integer :: i
      real,dimension(:) :: a
      real,dimension(size(a),size(a)) :: diagonal
      diagonal = 0
      do i = 1,size(a)
         diagonal(i,i) = a(i)
      end do
      ! return
end function diag


function find(a) result(found)
      integer i,j 
      logical,dimension(:) :: a
      integer,dimension(count(a)) :: found
      
      j = 1
      do i = 1,size(a)
        if (a(i)) then
            found(j) = i
            j = j + 1
        end if
      end do
      ! return
end function find

real function modularity(A)
    !    private :: B,k
        real::m,Q,delta_Q
        integer:: i,j,n_groups,n_current
        real,allocatable,dimension(:) :: s
        integer,dimension(:,:):: A
        integer,allocatable,dimension(:):: take,take_tmp,take2
        real,dimension(size(A,dim=1)):: k,g
        real,dimension(size(k),size(k)):: B
        real,allocatable,dimension(:,:):: Bg
        k = sum(A,dim=1)
        m = real(sum(k))/2.0
        modularity = 0
        !B = (/ (/ 1,2 /),(/ 1,2 /) /)
        g = 1
        
        
        do i = 1,size(k)
                do j = 1,size(k)
                        !print *,B(i,j)
                        B(i,j) = A(i,j)-(k(i)*k(j)/(2*m))
                end do
        end do
        
        n_groups = 1
        n_current = 1
        Q = 0
!        
        !endless loop
        do 
                !binary indexing. how is this best done in fortran? pack?
                if (allocated(take)) then 
                        deallocate(take)
                        deallocate(Bg)
                        deallocate(s)
                 endif
                        
                allocate(take(count(g==n_current)))
                allocate(s(count(g==n_current)))
                allocate(Bg(count(g==n_current),count(g==n_current)))
                take = find(g==n_current)
                !print*,g
                !print*,"bla"
                Bg = B(take,take) - diag(sum(B(take,take),2))
                 
                call dQ(Bg,k,m,delta_Q,s)
                
                if (delta_Q > 0.0001) then
!                        if (allocated(take_tmp) then
!                                deallocate(take_tmp)
!                        endif
!                        allocate(take_tmp(count(s==1)))
!                        take_tmp = take(find(s == 1))
                        !take2 = take(find(s == -1))
                        
               
   
                        n_groups = n_groups + 1
                        g(take(find(s == 1))) = n_groups
                        Q = Q + delta_Q
                else
                        !take = 1,size(A,1)
                        n_current = n_current + 1
                        if (n_current > n_groups) then 
                                exit
                        end if
                end if
                
        end do
       modularity = Q 
end function modularity
        
        subroutine dQ(B,k,m,delta_Q,s)
        
        
        integer :: i,j,ii,nrot,maxind
        real,dimension(:):: k
        
        real,dimension(:,:) :: B
        
        real,dimension(size(B,1),size(B,1)) :: u
        real :: m,delta_q,summe,tsumme,ttsumme
        real,dimension(size(B,1)) :: beta,s,used,t,ts,tts


        !this needs to return the eigenvalues in beta and the leading eigenvector in u
        call jacobi(B,beta,u,nrot)
        
        maxind = maxloc(beta,dim=1)
        s = sign(1.0,u(:,maxind))

        delta_Q = 0.0
        
        summe = 0 
        ts = s
        used = 0 !same size as s
        do i = 1,size(B,1)
        summe = summe + beta(i)*(dot_product(s,u(:,i)))**2
        end do
        
        tsumme = summe
        !fine tuning
        
        do i = 1,size(s,1)

            do j = 1,size(used)
            if (used(j)==0) then
              !  forall(j=1:size(used),used==0)
                tts = ts
                tts(j) = -1*tts(j)
          
                
                ttsumme = 0
!                print *, tts
!                print *,u(:,ii)
!                pause
                do ii = 1,size(B,1)
                ttsumme = ttsumme + beta(ii)*(dot_product(tts,u(:,ii)))**2
                end do

                if (ttsumme > tsumme) then
                    
                    ts = tts
                    tsumme = ttsumme
                    
                end if
             endif
             end do
            !end forall
            !used(find(s /= ts)) = 1
            where (s /= ts) used = 1
            
            s = ts
            if (summe >= tsumme) then
                exit
            end if
            summe = tsumme

      end do

        delta_Q = summe/(4*m)
        
        end subroutine dQ

real function wire(A,active)
     integer,dimension(:,:) :: A
     integer,dimension(:) :: active
     integer :: i,j,n_active
     integer,dimension(count(active>0)) :: choose
     real,dimension(count(active>0)-n_inp-1,2) :: pos 
     integer,dimension(count(active>0)-n_inp-1,count(active>0)-n_inp-1) :: Q
     integer,dimension(count(active>0)-n_inp-1,n_inp+1) :: B
     real,dimension (count(active>0)-n_inp-1,count(active>0)-n_inp-1) :: Q_inv
     real,dimension (n_inp+1,2) :: s
     !uses alpha = 1 throughout

     n_active =  count(active>0) 
     choose = find(active>0)
     Q = -A(choose(n_inp+1:n_active-1),choose(n_inp+1:n_active-1))
     
!     do j = 1,size(Q,1)
!        print *,Q(j,:),"\n"
!     end do
!     print *,"------\n"
     
     do i = 1,n_active-n_inp-1
         Q(i,i) = Q(i,i) - sum(Q(i,:))
     end do

!     do j = 1,size(Q,1)
!        print *,Q(j,:),"\n"
!     end do
!     print *,"------\n"

     do i = 1,n_active-n_inp-1
         B(i,1:n_inp) = A(choose(i+n_inp),1:n_inp)
         B(i,n_inp+1) = A(choose(i+n_inp),n_gates)
         Q(i,i) = Q(i,i) + sum(B(i,:))
     end do
     Q_inv = real(Q)
    
     s = 0
     do i = 1,n_inp
        s(i,1) = real(i-1)/(n_inp-1)
     end do
     s(n_inp+1,1) = 0.5
     s(n_inp+1,2) = 1.0


     pos(:,1) = matmul(B,s(:,1))
     pos(:,2) = matmul(B,s(:,2))
!     do j = 1,size(Q,1)
!        print *,Q(j,:),"\n"
!     end do
!        
!     print *,"------\n"
     call gaussj(Q_inv,pos)
!     print *,Q_inv
     wire = 0

     do i = n_inp+1,n_active-1
     do j = n_inp+1,n_active-1
     
     wire = wire + 0.5*A(choose(i),choose(j))*((pos(i-n_inp,1)-pos(j-n_inp,1))**2 + (pos(i-n_inp,2)-pos(j-n_inp,2))**2 )
     end do

     do j = 1,n_inp+1
     wire = wire + B(i-n_inp,j)*((pos(i-n_inp,1)-s(j,1))**2 + (pos(i-n_inp,2)-s(j,2))**2 )
     end do
     end do
end function wire

end program main
