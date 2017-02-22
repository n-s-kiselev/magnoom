Program FastFourierTransform
implicit none
double precision, parameter     :: pi  = 3.141592653589793238
double precision, parameter     :: tpi = 6.283185307179586477

character(len=200)              :: input_file_name
character(len=200)              :: output_file_name = "FFT_OUTPUT.csv"
character(len=200)              :: string
character(len=12 ), allocatable :: args(:)
integer                         :: num_args, iarg, fios
integer                         :: n, i, j, k, step

double precision                :: sx,sy,sz,et
double precision, allocatable   :: input(:)
double precision, allocatable   :: output(:)
double precision, allocatable   :: input_time(:)

num_args = command_argument_count()
if (num_args>0) then
    allocate(args(num_args))
    do iarg = 1, num_args
		call get_command_argument(iarg,args(iarg))
    enddo
    input_file_name = args(1)
	write(*,'(A20,A12)') "Trying to read from ",input_file_name
	if (file_exists(input_file_name)/=0) then
	  	write(*,'(A12,I2)') "Fatal error:",file_exists(input_file_name)
	else
		n = file_line_number(input_file_name)
		n = n-1 ! if file contains column title
		write(*,'(A24,I7)') "Number of lines in file: ", n
		i=2
		j=1
		do while (i<n)
			i=i*2
			j=j+1
		enddo
		n=i/2
		write(*,'(A10,I7,A3,I2,A7)') "I'll take ", n,"=2^",j-1 ," lines!" 
		allocate(input(2*n), output(2*n),input_time(n))
		call read_csv_file(input_file_name,input,n,5)!5 means values from the third column (5-mz, 3-mx)
		do i=1,2*n
			output(i)=input(i)
		enddo
		call four1(output,n,1)
		output_file_name = "FFT_"//trim(input_file_name)
		call write_csv_file(output_file_name,input,output,n)
		deallocate(input, output)
	endif
else
	write(*,'(A34)') "Fatal error: no input file :(        "
	write(*,'(A34)') "I'll do simple test culculation...   "
	n=8192
	allocate(input(2*n), output(2*n))
	call test_input(input,n)
	call test_input(output,n)
	call four1(output,n,1)
	write(*,'(A34)') "See file FFT_test_data.csv           "
	output_file_name = "FFT_test_data.csv"
	call write_csv_file(output_file_name,input,output,n)
!!!!!!!!DFT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call test_input(input,n)
	call test_input(output,n)
	call four0(output,n,1)
	write(*,'(A34)') "See file DFT_test_data.csv           "
	output_file_name = "DFT_test_data.csv"
	call write_csv_file(output_file_name,input,output,n)
	deallocate(input, output)
endif


CONTAINS

subroutine test_input(data,num) 
implicit none
integer            :: num
double precision   :: data(2*num)
integer            :: i,k
real               :: T1,T2
	T1 =1.0/100.0
	T2 =1.0/60.0
	do i=1, 2*num, 2
		k=(i-1)/2
	    data(i)  =  -exp(-0.01*tpi*k*T2)*cos(tpi*k*T1)!*sin(tpi*k*T2)
 	    !data(i)  =  -exp(-0.01*tpi*k*T2)*sin(tpi*k*T1)
	    data(i+1) = 0.0
	enddo
end subroutine test_input


function file_exists(filename) 
implicit none
integer           :: file_exists
character(len=12) :: filename
integer           :: ios !Input/Output Status
	open(100,file=filename, iostat=ios, action='read', status='old')
	close(100)
	file_exists = ios
end function file_exists



function file_line_number(filename) 
implicit none
integer            :: file_line_number
character(len=12)  :: filename
integer            :: ios,i
character(len=200) :: string
	open(100,file=filename, iostat=ios, action='read', status='old')
	i=0
	do while (ios/=-1)!while "not end of file"
	   read(unit=100,fmt='(A)', iostat=ios) string; i=i+1;
	enddo
	close(100)
	file_line_number = i
end function file_line_number



subroutine read_csv_file(filename,data,num,column) 
implicit none
integer            :: num
double precision   :: data(2*num)
integer            :: column
character(len=12)  :: filename
integer            :: ios,i,k
double precision   :: s(6)
character(len=200) :: string
	open(100,file=filename, iostat=ios, action='read', status='old')
	read(100,*) string
	i=0
	do i=1, 2*num, 2
	   read(100,*) s(1),s(2),s(3),s(4),s(5),s(6)
	   data( i ) = s(column) ! real part
	   data(i+1) = 0.0       ! imaginary part
	   input_time((i+1)/2) = s(2)!time column
	enddo
	close(100)
end subroutine read_csv_file



subroutine write_csv_file(filename,datain,dataout,num) 
implicit none
integer            :: num
double precision   :: datain(2*num)
double precision   :: dataout(2*num)
character(len=200) :: filename
integer            :: ios,i,k
double precision   :: omega
	open(100,file=filename, iostat=ios, action='write', status='unknown')
 	!write(100,'(A17)') "n, signal, Im, Re"
 	write(100,'(A29)') " time, signal, omega, Im, Re,"
	i=0
	do i=1, 2*num, 2
		k = (i-1)/2
	    !write(100,'(I6,5(A1,ES14.7))') &
	    !k,",", input(i),",", output(i)*k,",", output(i+1)*k,",",sqrt(output(i)**2),",", sqrt(output(i+1)**2)
	    omega = 6.283185307179586*k/(input_time(num))
	    write(100,'(5(ES14.7,A1))') &
	    input_time(k+1),",", input(i),",",omega,",", output(i)*k,",", output(i+1)*k,","
	enddo
	close(100)
end subroutine write_csv_file


!The following subroutine is taken from "Numerical Recipes In F77" ISBN 0-521-43064-X

subroutine four1(data,nn,isign)
implicit none
integer          :: isign,nn
double precision :: data(2*nn)
integer          :: i,istep,j,m,mmax,n
real             :: tempi,tempr
double precision :: theta,wi,wpi,wpr,wr,wtemp
	n=2*nn
	j=1
	do i=1,n,2
		if(j.gt.i)then
			tempr=data(j)
			tempi=data(j+1)
			data(j)=data(i)
			data(j+1)=data(i+1)
			data(i)=tempr
			data(i+1)=tempi
		endif
		m=n/2
1 		if ((m.ge.2).and.(j.gt.m)) then
			j=j-m
			m=m/2
			goto 1
		endif
		j=j+m
	enddo
	mmax=2
2 	if (n.gt.mmax) then 
		istep=2*mmax
		theta=6.28318530717959d0/(isign*mmax)
		wpr=-2.d0*sin(0.5d0*theta)**2
		wpi=sin(theta)
		wr=1.d0
		wi=0.d0
		do  m=1,mmax,2
			do  i=m,n,istep
				j=i+mmax
				tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
				tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
				data(j)=data(i)-tempr
				data(j+1)=data(i+1)-tempi
				data(i)=data(i)+tempr
				data(i+1)=data(i+1)+tempi
			enddo 
			wtemp=wr
			wr=wr*wpr-wi*wpi+wr
			wi=wi*wpr+wtemp*wpi+wi
		enddo 
		mmax=istep
		goto 2
	endif
	return
end subroutine four1

subroutine four0(data_in,nn,isign)
implicit none
integer          :: isign,nn
double precision :: data_in(2*nn)
double precision :: data_outR(nn)
double precision :: data_outI(nn)
integer          :: i,k,ii,kk,n
	n=2*nn
	do i=1,n,2
	ii=(i-1)/2
	data_outR(ii+1)=0
	data_outI(ii+1)=0
		do k=1,n,2
		kk=(k-1)/2
		data_outR(ii+1) = data_outR(ii+1)+data_in(k)*cos(6.28318530717959d0*ii*kk/nn)
		data_outI(ii+1) = data_outI(ii+1)+data_in(k)*sin(6.28318530717959d0*ii*kk/nn)
		enddo	
	enddo
	do i=1,n,2
		ii=(i-1)/2
		data_in(i)=data_outR(ii+1)
		data_in(i+1)=data_outI(ii+1)
	enddo
	return
end subroutine four0

end program FastFourierTransform


