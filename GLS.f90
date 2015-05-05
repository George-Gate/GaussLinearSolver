!GLS.f90 
!
! Windows下若要开启并行功能，请加上编译选项 /Qopenmp
! 或者在VS中进行如下设置：项目属性->Fortran->Language->Process OpenMP Directives: Generate Parallel Code (/Qopenmp)
! 建议不要使用Intel的超线程技术，以下是某次运行的测试结果(CPU为Core i7 4702MQ, 4核心8线程)
! 串行：     2:55.03
! 并行4线程：2:16.99
! 并行8线程：2:34.27
! 此程序目前的瓶颈在内存访问上，并行加速比最大是1.3左右，待优化。
!
![Modules]
!  GLS - 用高斯消元法解线性方程组Ax=b
!  
![Functions]
!  GaussLinearSolver - 用高斯消元法解线性方程组Ax=b
!       [Type 1]
!       function GaussLinearSolver(A,b,n) result(x)
!       Input: A - 系数矩阵，一个n阶可逆矩阵，A(i,j)表示第j行第i列的元素
!              b - 常数向量，一个n维列向量
!              n - 方程组未知数的个数
!       Result: x - 方程组Ax=b的解
!       Memory Usage: Size(A)+3*Size(b)。 相较于[Type 2]形式的调用，这种调用方式会耗费更多的内存。
!       
!       [Type 2]
!       function GaussLinearSolver(Ab,n) result(x)
!       Input: Ab - 方程组的增广矩阵，要求R(Ab)=n，Ab(i,j)表示第j行第i列的元素
!              n - 方程组未知数的个数，Ab的尺寸为A(n+1,n)
!       Result: x - 方程组Ax=b的解，一个n维列向量
!       Memory Usage: Size(b)，在输入数组上就地消元。
!

    module GLS
        private

        public GaussLinearSolver

        interface GaussLinearSolver
            module procedure GaussLinearSolver_A
            module procedure GaussLinearSolver_Ab
        end interface GaussLinearSolver
        
    contains
        function GaussLinearSolver_A(A,b,n) result(x)
            implicit none
            integer,intent(in) :: n
            real(kind=8),intent(in)  :: A(n,n),b(n)
            real(kind=8),allocatable :: x(:)
            real(kind=8),allocatable :: Ab(:,:)
            
            allocate(x(n))
            allocate(Ab(n+1,n))
            Ab(1:n,:)=A(:,:)
            Ab(n+1,:)=b(:)
            x=GaussLinearSolver_Ab(Ab,n)
            
            return
        end function GaussLinearSolver_A
        
        function GaussLinearSolver_Ab(Ab,n) result(x)
        !用高斯消元法解线性方程组Ax=b
        !Ab是方程组的增广矩阵,n是未知数个数
        !Ab(i,j)表示第j行第i列的元素
            use ifport
            implicit none
            !设置计算精度
            real(kind=8),parameter :: eps=1.0d-4
            integer, parameter :: cmdWidth=80
            
            integer,intent(in) :: n
            real(kind=8),intent(inout) :: Ab(n+1,n)
            real(kind=8),allocatable :: x(:)
            integer :: i,j,maxJ, progressCounter
            real(kind=8) :: maxOfCol
            real(kind=8) :: time1

            time1=dclock()       ! 计时
            write (unit=6,fmt="(A,I0)") "@GaussLinearSolver:  N=",n
            write (unit=6,fmt="(A)") "Progress:"
            write (unit=6,fmt="(<cmdWidth-1>A)") ("*",i=1,cmdWidth-1)
            write (unit=6,fmt="(<cmdWidth-7>A)") "0%",(" ",i=1,(cmdWidth-10)/2+mod(cmdWidth,2)),"50%",(" ",i=1,(cmdWidth-10)/2),"100%"
            progressCounter=0
            
            allocate(x(n))
            
            do j=1,n
                maxOfCol=0d0
                do i=j,n                        !从第j行开始往下找第j列的最大值
                    if (abs(Ab(j,i))>maxOfCol) then
                        maxOfCol=abs(Ab(j,i))
                        maxJ=i
                    end if
                end do

                if (maxOfCol<eps) then          !若最大值为零，则无法确保主对角元素非零，矩阵A不可逆
                    write (unit=6,fmt="(A)") "@GaussLinearSolver[Function]: error:输入的A不是可逆矩阵！"
                    x=0d0
                    return
                else if (maxJ/=j) then
                    if (Ab(j,j)*Ab(j,maxJ)>=0) then   !将maxJ行倍加至第j行, 先判断主元符号是否相同
                        !dir$ parallel
                        Ab(j:n+1,j)=Ab(j:n+1,j)+Ab(j:n+1,maxJ)
                    else
                        !dir$ parallel
                        Ab(j:n+1,j)=Ab(j:n+1,j)-Ab(j:n+1,maxJ)
                    end if
                end if
                
                !dir$ parallel
                Ab(j:n+1,j)=Ab(j:n+1,j)/Ab(j,j)              !单位化
                

                !$omp parallel default(none) shared(Ab,j,n)
                !$omp do schedule(guided) private(i)
                    do i=j+1,n                                !消元：用第j行将A(j,j)下方的元素归零,最终得到上三角矩阵
                        !dir$ parallel
                        Ab(j+1:n+1,i)=Ab(j+1:n+1,i)-Ab(j+1:n+1,j)*Ab(j,i)
                        Ab(j,i)=0d0                           !防止出现1.0E-38这样的数......
                    end do
                !$omp end do
                !$omp end parallel

                ! 输出进度
                do while (progressCounter/real(cmdWidth-1)<(real(j)*(j**2d0-3d0*j*n+3d0*n**2d0-1d0))/(n*(n-1d0)*(n+1d0)))
                    i=putc('*')
                    progressCounter=progressCounter+1
                end do
            end do

            !dir$ parallel
            x(:)=Ab(n+1,:)
            
            do j=n,2,-1               ! 将上三角矩阵消元成对角阵,使用的对角元为A(j,j)
                do i=j-1,1,-1
                    x(i)=x(i)-Ab(j,i)*x(j)
                    Ab(j,i)=0d0
                end do
            end do
            
            Ab(n+1,:)=x(:)
            
            write (unit=6,fmt="(<cmdWidth-1-progressCounter>A)") ("*",i=progressCounter+1,cmdWidth-1)
            write (unit=6,fmt="(A,F0.3,A)") "Solving procedure finished.(Elapsed time:",dclock()-time1,"s)"
        
            return
        end function GaussLinearSolver_Ab
    end module GLS

    