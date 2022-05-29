module gen_dtree
  use var_global
  implicit none

  public :: Node

  type Node
     integer depth
     integer pclass 
     integer :: findex = 0
     real :: threshold = 0.0
     real :: gini = 1.0
     !type(Node), pointer :: parent
     type(Node), pointer :: left => NULL()
     type(Node), pointer :: right => NULL()
  end type Node


contains

!  subroutine tree_1D(set1, set2, numb, area)
!    real*8 set1(:), set2(:), area
!    real(8), dimension(size(set1)+size(set2), 1) :: X
!    integer, dimension(size(set1)+size(set2)) :: y
!    integer i, numb, ns1, n
!
!    ns1 = size(set1)
!    n = (size(set1)+size(set2))
!
!    do i=1, n
!       if (i .le. ns1) then
!          X(i, 1) = set1(i)
!          y(i) = 1
!       else
!          X(i, 1) = set2(i-ns1)
!          y(i) = 2
!       end if
!    end do
!
!    call fit_tree(X, y, max_depth_1D, area, numb)
!
!  end subroutine tree_1D
!  


!  subroutine tree_2D(set1, set2, width, numb, area)
!    real*8 set1(:,:), set2(:,:), area, width
!    real(8), dimension((size(set1)+size(set2))/2, 2) :: X
!    integer numb, i, ns1, n
!    integer, dimension((size(set1)+size(set2))/2) :: y
!
!    ns1 = size(set1)/2
!    n = (size(set1)+size(set2))/2
!
!    do i=1, n
!       if (i .le. ns1) then
!          X(i, :) = set1(i, :)
!          y(i) = 1
!       else
!          X(i, :) = set2(i-ns1, :)
!          y(i) = 2
!       end if
!    end do
!
!    call fit_tree(X, y, max_depth_2D, area, numb)
!    width = area
!    
!  end subroutine tree_2D
!

  subroutine di_tree_1D(X, ngroups, overlap_n, cvacc)
    real*8 X(:), cvacc
    integer max_depth, cvmc, ngroups(:), i, cl , j, overlap_n
    integer, dimension(size(X)) :: y
    real(8), dimension(size(X), 1) :: Xn

    cvacc = 0.0
    cvmc = 0
    cl = 1

    do i=1, 999
       if (ngroups(i) .ne. 0) then

          if (i .eq. 1) then
             do j=1, ngroups(i)
                y(j) = cl
             end do
          else
             do j=sum(ngroups(:i-1))+1, sum(ngroups(:i))
                y(j) = cl
             end do
          end if
          
          cl = cl + 1
       end if
    end do

    do i=1, size(X)
       Xn(i, :) = X(i)
    end do

    call fit_tree(Xn, y, max_depth_1D, CV_fold, cvacc, cvmc)
    
    overlap_n = cvmc

  end subroutine di_tree_1D


  subroutine sis_tree_1D(X, ngroups, overlap_n)
    real*8 X(:), cvacc
    integer max_depth, cvmc, ngroups(:), i, cl , j, overlap_n
    integer, dimension(size(X)) :: y
    real(8), dimension(size(X), 1) :: Xn
    logical feat_invalid

    feat_invalid = .false.

    cvacc = 0.0
    cvmc = 0
    cl = 1

    do i=1, 999
       if (ngroups(i) .ne. 0) then

          if (i .eq. 1) then
             do j=1, ngroups(i)
                y(j) = cl
             end do
          else
             do j=sum(ngroups(:i-1))+1, sum(ngroups(:i))
                y(j) = cl
             end do
          end if

          cl = cl + 1
       end if
    end do

    do i=1, size(X)
       Xn(i, :) = X(i)
    end do

    ! check for NaN                                                                                                     
    do i=1, size(X)
       if (isnan(X(i)) .or. (X(i) .lt. maxfval_lb) .or. (X(i) .gt. maxfval_ub) .or. all_same(X)) then
          feat_invalid = .true.
       end if
    end do

    if (feat_invalid) then
       overlap_n = size(y)
    else
       call fit_tree(Xn, y, max_depth_1D, 1, cvacc, cvmc)
       overlap_n = cvmc
    end if
    
  end subroutine sis_tree_1D  


  subroutine di_tree_2D(X, ngroups, overlap_n, cvacc)
    real*8 X(:,:), cvacc
    integer max_depth, cvmc, ngroups(:), i, cl , j, overlap_n
    integer, dimension(int(size(X)/2)) :: y

    cvacc = 0.0
    cvmc = 0
    cl = 1

    do i=1, 999
       if (ngroups(i) .ne. 0) then

          if (i .eq. 1) then
             do j=1, ngroups(i)
                y(j) = cl
             end do
          else
             do j=sum(ngroups(:i-1))+1, sum(ngroups(:i))
                y(j) = cl
             end do
          end if
          
          cl = cl + 1
       end if
    end do


    call fit_tree(X, y, max_depth_2D, CV_fold, cvacc, cvmc)
    
    overlap_n = cvmc

  end subroutine di_tree_2D



  subroutine fit_tree(X, y, max_depth, cv, cvacc, cvmc)
    real*8 X(:,:), acc, cvacc
    integer y(:), max_depth, nfeatures, n, c, i, k, numb, depth, p, nsamples, ncorrect, cv, fs, ista, iend, cvmc
    type(Node), pointer :: pnode, rnode
    integer , dimension(2) :: xshape
    integer, dimension(size(y)) :: rand_inds
    real(8), dimension(:,:), allocatable :: xtrain, xval
    integer, dimension(:), allocatable :: train_inds, val_inds, yval, ytrain

    nsamples = size(y)
    fs = int(float(nsamples)/float(cv))
    rand_inds = randomized_indices(size(y))
    xshape = shape(X)
    nfeatures = xshape(2)
    cvmc = 0
    cvacc = 0.0

    do k=1, cv
       ista = max(1, int((k-1)*fs)+1)

       if (cv .eq. 1) then
          ista = 1
       end if
       
       iend = min(nsamples, (k*fs))
       
       allocate(val_inds(iend-ista+1))
       allocate(yval(iend-ista+1))
       allocate(xval(iend-ista+1, nfeatures))

       if (cv .gt. 1) then
          allocate(train_inds(nsamples-(iend-ista+1)))
          allocate(ytrain(nsamples-(iend-ista+1)))
          allocate(xtrain(nsamples-(iend-ista+1), nfeatures))
       else
          allocate(train_inds(nsamples))
          allocate(ytrain(nsamples))
          allocate(xtrain(nsamples, nfeatures))
       end if
       

       c = 1
       do i=ista, iend
          val_inds(c) = rand_inds(i)
          c = c + 1
       end do

       if (cv .gt. 1) then
          c = 1
          do i=1, nsamples
             if (.not. any(val_inds .eq. i)) then
                train_inds(c) = i
                c = c+ 1
             end if
          end do
       else
          c = 1
          do i=1, nsamples
             train_inds(c) = i
             c = c+ 1
          end do
       end if
       
       

       do i=1, size(val_inds)
          xval(i, :) = X(val_inds(i), :)
          yval(i) = y(val_inds(i))
       end do

       do i=1, size(train_inds)
          xtrain(i, :) = X(train_inds(i), :)
          ytrain(i) = y(train_inds(i))
       end do


       depth = 1


       n = 0
       c = -1

       do i=1, size(y)
          if (y(i) .ne. c) then
             n = n + 1
             c = y(i)
          end if
       end do


       call grow_tree(xtrain, ytrain, nfeatures, n, pnode, depth, max_depth)
       allocate(rnode)
       rnode = pnode


       ncorrect = 0

       do i=1, size(yval)

          p = 0
          call pred_one(Xval(i, :), p, pnode)
          pnode = rnode
          if (p .eq. yval(i)) then
             ncorrect = ncorrect + 1
          end if
       end do

       acc = float(ncorrect)/float(size(yval))

       numb = size(yval) - ncorrect

       cvmc = cvmc + numb
       cvacc = cvacc + acc

       call free_tree(pnode)
       deallocate (rnode)

       deallocate(val_inds)
       deallocate(train_inds)
       deallocate(yval)
       deallocate(ytrain)
       deallocate(xval)
       deallocate(xtrain)

    end do

    !cvmc = int(float(cvmc)/float(cv))
    cvacc = cvacc/float(cv)
    cvacc = 1.0 - cvacc

  end subroutine fit_tree


  recursive subroutine free_tree(tree)
    type(Node), pointer :: tree

    if (.not. associated(tree)) then
       return
    end if

    if (associated(tree%left)) then
       call free_tree(tree%left)
    end if
    
    if (associated(tree%right)) then
       call free_tree(tree%right)
    end if

    deallocate(tree)
    nullify(tree)
    
  end subroutine free_tree


  recursive subroutine pred_one(X, p, cnode)
    real*8 x(:)
    integer p
    type(Node), pointer :: cnode

    if (.not. associated(cnode)) then
       return
    end if

    p = cnode%pclass

    !leaf
    if (cnode%findex .eq. -1) then
       return
    end if
    
    if (x(cnode%findex) .le. cnode%threshold) then
       call pred_one(X, p, cnode%left)
    else
       call pred_one(X, p, cnode%right)
    end if

  end subroutine pred_one
  
  

  recursive subroutine predict_one(X, p, cnode)
    real*8 x(:)
    integer p, c, x2
    type(Node), pointer:: cnode
    type(Node), pointer :: temp

    x2 = size(X)
    print *, p, X, cnode%depth, cnode%findex, cnode%threshold
    if (associated(cnode%left) .and. (cnode%findex .gt. 0) ) then
       if ((X(cnode%findex) .lt. cnode%threshold) .and. (cnode%pclass .ne. -1)) then

          cnode = cnode%left
          if (cnode%pclass .ne. -1) then
             p = cnode%pclass
          end if

          call predict_one(X, p, cnode)

       end if
    end if
    
    if (associated(cnode%right) .and. (cnode%findex .gt. 0) .and. (cnode%pclass .ne. -1)) then
       cnode = cnode%right

       if (cnode%pclass .ne. -1) then
          p = cnode%pclass
       end if

       call predict_one(X, p, cnode)

    end if

  end subroutine predict_one



  recursive subroutine grow_tree(X, y, nfeatures, nclasses, cnode, depth, mdepth)

    type(Node), pointer :: cnode
    integer :: depth
    real*8 X(:, :), gini
    integer y(:), nfeatures, nclasses, ymin, mdepth
    
    integer, dimension(size(y)) :: ind_next
    
    real(8), allocatable, dimension(:,:) :: xleft, xright
    integer, allocatable, dimension(:) :: yleft, yright
    integer, dimension(nclasses) :: classes

    integer pclass, ibest_feature, i, j, nleft, nright, cmax
    real(8) loss, threshold

    threshold = 0.0
    loss = 1.0
    ibest_feature = -1
    cmax = 0
    gini = 0.0

    !adjust y vals if required
    ymin = minval(y) - 1

    do i=1, nclasses
       classes(i) = 0
    end do

    do i=1, size(y)
       classes(y(i)) = classes(y(i)) + 1
    end do

    do i=1, nclasses
       if (classes(i) .ge. cmax) then
          cmax = classes(i)
          pclass = i
       end if
    end do

    do i=1, nclasses
       gini = gini + ((float(classes(i))/float(size(y)))**2.0)
    end do
    gini = 1.0 - gini


    allocate(cnode)
    cnode = Node(depth=depth, pclass=pclass)
    if ((gini .gt. 0.001) .and. (depth .le. mdepth)) then

       call gini_loss(X, y, nfeatures, nclasses, threshold, loss, ibest_feature)


       cnode%findex = ibest_feature
       cnode%threshold = threshold
       cnode%gini = loss

       nleft = 0
       nright = 0

       do i=1, size(y)
          if (X(i, ibest_feature) .le. threshold) then
             nleft = nleft + 1
             ind_next(i) = 1
          else
             nright = nright + 1
             ind_next(i) = 2
          end if
       end do

       if (nleft .gt. 2) then
          allocate(xleft(nleft, nfeatures))
          allocate(yleft(nleft))

          j = 1
          do i=1, size(y)
             if (ind_next(i) .eq. 1) then
                xleft(j, :) = X(i, :)
                yleft(j) = y(i)
                j = j + 1
             end if
          end do

          call grow_tree(xleft, yleft, nfeatures, nclasses, cnode%left, depth+1, mdepth)

          deallocate(xleft)
          deallocate(yleft)
       end if
          

       if (nright .gt. 2) then
          allocate(xright(nright, nfeatures))
          allocate(yright(nright))

          j = 1
          do i=1, size(y)
             if (ind_next(i) .eq. 2) then
                xright(j, :) = X(i, :)
                yright(j) = y(i)
                j = j + 1
             end if
          end do

          call grow_tree(xright, yright, nfeatures, nclasses, cnode%right, depth+1, mdepth)

          deallocate(xright)
          deallocate(yright)
       end if

    end if

  end subroutine grow_tree


  function all_same(X) result(same)
    real*8 X(:)
    logical :: same

    same = .true.

    if (sum(abs(X - X(1))) .gt. 1d-9) then
       same = .false.
    end if

  end function all_same
  
  subroutine sis_gini_1D(X, ngroups, overlap_n)
    real*8 X(:), gini, thr
    integer max_depth, ngroups(:), i, cl , j, overlap_n, ibest
    integer, dimension(size(X)) :: y
    real(8), dimension(size(X), 1) :: Xn
    logical feat_invalid

    gini = 1.0
    feat_invalid = .false.

    cl = 1

    do i=1, 999
       if (ngroups(i) .ne. 0) then

          if (i .eq. 1) then
             do j=1, ngroups(i)
                y(j) = cl
             end do
          else
             cl = cl + 1
             do j=sum(ngroups(:i-1))+1, sum(ngroups(:i))
                y(j) = cl
             end do
          end if
       end if
    end do

    do i=1, size(X)
       Xn(i, :) = X(i)
    end do

    ! check for NaN
    do i=1, size(X)
       if (isnan(X(i)) .or. (X(i) .lt. maxfval_lb) .or. (X(i) .gt. maxfval_ub) .or. all_same(X)) then
          feat_invalid = .true.
       end if
    end do
    
    if (feat_invalid) then
       overlap_n = size(y)
    else

       call gini_loss(Xn, y, 1, cl, thr, gini, ibest)
       overlap_n = int(gini*float(size(y)))

    end if
    
  end subroutine sis_gini_1D

  
  subroutine gini_loss(X, y, nfeatures, nclasses, bf_threshold, loss, ibest_feature)
    ! X is (N x M) feature matrix
    ! y is labels - can be more than 2

    real*8 X(:,:), threshold, bf_threshold, gleft, gright, gini, gini_best, gini_best_all, loss
    integer i, j, k, c, nsamples, nclasses, nfeatures, ibest_feature, y(:), ymin
    real(8), dimension(size(y)) :: xc, xs
    integer, dimension(nclasses) :: cleft, cright
    integer, dimension(size(y)) :: ys, indsarray

    nsamples = size(y)
    ibest_feature = -1
    gini_best_all = 1.0
    bf_threshold = -1.0
    ymin = minval(y) - 1

    
    do i=1, nfeatures

       gini_best = 1.0
       gini = 1.0

       xc = X(:, i)
       indsarray = argsort(xc)

       do j=1, nclasses
          cleft(j) = 0
          cright(j) = 0
       end do

       do j=1, nsamples
          cright(y(j)-ymin) = cright(y(j)-ymin) + 1
       end do

       do j=1, nsamples
          xs(j) = xc(indsarray(j))
          ys(j) = y(indsarray(j))
       end do

       do j=2, nsamples
          c = ys(j-1)
          cleft(c-ymin) = cleft(c-ymin) + 1
          cright(c-ymin) = cright(c-ymin) - 1

          gleft = 0.0
          gright = 0.0

          do k=1, nclasses
             gleft = gleft + ( (float(cleft(k)) / float(j))**2.0 )
             gright = gright + ( (float(cright(k)) / float(nsamples-j+1))**2.0 )
          end do

          gleft = 1.0 - gleft
          gright = 1.0 - gright

          if ((nsamples-j .ne. 0) .and. (gright .ne. gright)) then
             write(*,*) "GRIGHT NAN", cleft, cright, nsamples, j
          end if

          gini = ((gleft * float(j)) + (gright * float(nsamples-j+1))) / float(nsamples)

          if (gini .lt. 0.000000001) then
             write(*,*) "NEG GINI", cleft, cright, nsamples, j
          end if

          if (gini .lt. gini_best) then
             gini_best = gini
             if (xs(j) .ne. (xs(j-1))) then
                threshold = (xs(j) + xs(j-1)) / 2.0
             end if
          end if
       end do

       if (gini_best .lt. gini_best_all) then
          gini_best_all = gini_best
          ibest_feature = i
          bf_threshold = threshold
       end if
    end do
    
    loss = gini_best_all

  end subroutine gini_loss


  function randomized_indices(m) result(ind)
    integer m, i, j, k, temp
    integer, dimension(m) :: ind

    ind = (/(i, i=1, m)/)
    do i=1, m-1
       j = m - i
       k = int((rand(0) * float(j)) + 1)
       temp = ind(j)
       ind(j) = ind(k)
       ind(k) = temp
    end do
  end function randomized_indices

  function argsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    real(8), intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    real(8) :: temp2
    real(8) :: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
       ! find ith smallest in 'a'
       imin = minloc(a2(i:),1) + i - 1
       ! swap to position i in 'a' and 'b', if not already there
       if (imin /= i) then
          temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
           temp1 = b(i); b(i) = b(imin); b(imin) = temp1
       end if
     end do
   end function argsort
  
end module gen_dtree
