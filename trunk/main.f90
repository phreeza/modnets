program main

use nr,only:indexx
implicit none

!integer,parameter :: n_pop=1000,L=300,epoch=100,n_gen=100000,n_gates=30,n_gates_nom=20,n_transient=5,n_inp=4
integer,parameter :: n_pop=1000,L=300,epoch=100,n_gen=100000,n_gates=50,n_gates_nom=30,n_transient=5,n_inp=8
integer :: goal,gen,gate1,gate2
integer :: a,b,c
integer,dimension(n_inp , 2**n_inp) :: inp = 0
integer,dimension(n_gates) :: active = 0
integer,dimension(n_gates,n_pop) :: x=0,y=0
integer,dimension( 2**n_inp) :: out = 0
integer,dimension( n_pop ) :: order = 0
integer,dimension(n_gates , n_pop):: in1=0 , in2=0
integer,dimension(n_gates , n_pop) :: state=0
real, dimension(n_pop) :: qual = 0.0
real,dimension(4) :: rand
logical :: again = .true.
integer old_sum
CALL RANDOM_SEED
!CALL RANDOM_SEED 
!CALL RANDOM_SEED(put = (/ 123, 222 /))
!call random_number(rand)

do a = 1,2**n_inp
   do b = 1,n_inp
      inp(b,a)=transfer(btest(a-1,b),0)
   end do
end do

out = iand(ior(ieor(inp(1,:),inp(2,:)),ieor(inp(3,:),inp(4,:))),ior(ieor(inp(5,:),inp(6,:)),ieor(inp(7,:),inp(8,:))))
!out = iand(ieor(inp(1,:),inp(2,:)),ieor(inp(3,:),inp(4,:)))
!out(5) = 0
print*, out 
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
      !if (sum(active(n_inp+1:n_gates-1)) /= n_gates_nom) then
          !   again = .true.
      !else
         again = .false.
          !end if
   end if 

end do
!end net generation loop
end do

print*,"Go!"

!start main loop
do gen = 1,n_gen


print*,gen
qual = 0

!evaluate nets
do c = 1,2**n_inp
    state = 0
        x = 0
        y = 0
        do a = 1,n_transient
                state(1:n_inp,:) = spread(inp(:,c),ncopies=n_pop,dim=2)
                do b = 1,n_pop
                        x(n_inp+1:n_gates,b) = state(in1(n_inp+1:n_gates,b),b)
                        y(n_inp+1:n_gates,b) = state(in2(n_inp+1:n_gates,b),b)
                end do
                state(n_inp+1:n_gates,:) = 1 - x(n_inp+1:n_gates,:)*y(n_inp+1:n_gates,:)
        end do
        

        qual = qual +  NOT(IEOR(state(n_gates,:),spread(out(c),ncopies=n_pop,dim=1)))+2
                
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
          
   !else
      !if (sum(active(n_inp+1:n_gates-1)) /= n_gates_nom) then
          !   qual(a) = 0.0
          !end if
   end if 
end do



!sort nets by quality

call indexx(qual,order)
print*,maxval(qual)/2.**n_inp,real(sum(qual))/(2.**n_inp*n_pop)
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

!contains 
!
!function modularity(A)
!        k = sum(A)
!        m = real(sum(k))
!        B = (/ (/ 1,2 /),(/ 1,2 /) /)
!        g = 1
!
!        
!        do i = 1,size(k)
!                do j = 1,size(k)
!                        B(i,j) = A(i,j)-(k(i)*k(j)/(2*m))
!                end do
!        end do
!        
!        n_groups = 1
!        n_current = 1
!        Q = 0
!        
!        !endless loop
!        do
!        
!        
!                !binary indexing. how is this best done in fortran? pack?
!                take = find(g==n_current)
!                Bg = B(take,take) - diag(sum(B(take,take),2))
!                
!                call dQ(Bg,k,m,delta_Q,s)
!                
!                if (delta_Q > 0.0001) then
!                        take_tmp = take(find(s == 1))
!                        take2 = take(find(s == -1))
!                        take = take_tmp
!                        n_groups = n_groups + 1
!                        Q = Q + delta_Q
!                else
!                        take = 1,size(A,1)
!                        n_current = n_current + 1
!                        if (n_current > n_groups) then exit
!                end if
!                
!        end do
!        
!        
!        
!        contains subroutine dQ(B,k,m,deltaQ,s)
!        !this needs to return the eigenvalues in beta and the leading eigenvector in u
!        call eigen(B,u,beta)
!        s = sign(u)
!        
!        ts = s
!        used = 0 !same size as s
!        do i = 1,size(B,1)
!            summe = summe + beta(i)*(s*u(:,i))^2
!        end do
!        
!        tsumme = summe
!        !fine tuning
!        
!        do i = 1,size(s,1)
!
!      !      do j = find(used == 0)
!                forall(j=1:size(used),used==0)
!                tts = ts
!                tts(j) = -1*tts(j)
!          
!                
!                ttsumme = 0
!
!                do i = 1,size(B,1)
!                    ttsumme = ttsumme + beta(i)*(tts*u(:,i))^2
!                end do
!
!                if (ttsumme > tsumme) then
!                    
!                    ts = tts
!                    tsumme = ttsumme
!                    
!                end if
!            end forall
!            !used(find(s /= ts)) = 1
!            where (s /= ts) used = 1
!            
!            s = ts
!            if (summe >= tsumme) then
!                break
!            end if
!            summe = tsumme
!
!      end do
!
!        delta_Q = summe/(4*m)
!        
!        end subroutine dQ
!
!end function modularity


end program main
