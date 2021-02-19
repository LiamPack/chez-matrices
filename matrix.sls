(library (matrix) (export make-matrix matrix-ref matrix-set! matrix-dims matrix-num-dims
          matrix-num-vals matrix-rows matrix-cols matrix? matrix-copy
          matrix-flatten vec->matrix matrix-ref-row matrix-set-row!
          matrix-set-diagonal! matrix-identity
          matrix-max matrix-min matrix-contract matrix-ref-col matrix-set-col!
          matrix-map mul T matrix-fold-left tr euclidean-norm cross-product
          theta phi linsolve matrix-inverse matrix= make-matrix do-matrix matrix-fold)
  (import (rnrs))

  (define (matrix-dims tensor)
    (let tensor-dims-rec ([tensor tensor] [dim-list (list)])
      (if (number? (vector-ref tensor 0))
          (reverse (cons (vector-length tensor) dim-list))
          (tensor-dims-rec (vector-ref tensor 0) (cons (vector-length tensor) dim-list)))))

  (define (matrix-num-dims tensor)
    (length (matrix-dims tensor)))

  (define-syntax do-matrix
    (syntax-rules ()
      [(_ m () rest ...) (begin rest ...)]
      [(_ m (() j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i 1))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m (i j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i (car (matrix-dims m))))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m ret (i j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i (car (matrix-dims m))) ret)
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]))

  (define-syntax make-matrix
    (syntax-rules ()
      [(_ i) (make-vector i)]
      [(_ i j)
       (do ([m (make-vector i)]
            [k 0 (+ k 1)])
           ((= k i) m)
         (vector-set! m k (make-vector j)))]
      [(_ i j k ...)
       (make-vector i (make-matrix j k ...))]))

  (define-syntax matrix-ref
    (syntax-rules ()
      [(_ m i) (vector-ref m i)]
      [(_ m i j k ...) (vector-ref (matrix-ref m i k ...) j)]))

  (define-syntax matrix-set!
    (syntax-rules ()
      [(_ m i val) (vector-set! m i val)]
      [(_ m i j k ... val) (matrix-set! (vector-ref m i) j k ... val)]))


  (define (matrix-num-vals m) (map * (matrix-dims m)))
  (define (matrix-rows a) (vector-length a))
  (define (matrix-cols a)
    (if (vector? (matrix-ref a 0))
        (vector-length (matrix-ref a 0))
        1))

  (define matrix?
    (lambda (x)
      (and (vector? x)
           (> (vector-length x) 0)
           (vector? (vector-ref x 0)))))

  ;; assert proper dims
  (define (matrix-copy a)
    (let* ([m (make-matrix (matrix-rows a) (matrix-cols a))])
      (do-matrix m m (i j) (matrix-set! m i j (matrix-ref a i j)))))

  (define (matrix-copy! a b)
    (do-matrix a (i j) (matrix-set! a i j (matrix-ref b i j))))

  (define (vec->matrix v)
    (make-matrix (vector-length v) 1))

  (define (matrix-flatten m)
    (let* ([v (make-vector (* (matrix-cols m) (matrix-rows m)))])
      (do-matrix m v (i j)
                 (vector-set! v (+ j (* i (matrix-cols m)))
                              (matrix-ref m i j)))))

  (define (matrix-ref-row m i) (vector-ref m i))
  (define (matrix-set-row! m i v) (vector-set! m i v))

  ;; verry inefficient but clean!
  (define (matrix-ref-col m i)
    (matrix-ref-row (T m) i))
  (define (matrix-set-col! m i v)
    (do-matrix m (row ()) (matrix-set! m row i (vector-ref v row))))

  ;; assert len(v) = min(dimensions(m))
  (define (matrix-set-diagonal! m v)
    (let ([new-m (matrix-copy m)])
      (do ([i 0 (+ 1 i)])
          ((>= i (vector-length v)) m)
        (matrix-set! m i i (vector-ref v i)))))

  (define (matrix-identity dim)
    (let ([m (make-matrix dim dim)])
      (matrix-set-diagonal! m (make-vector dim 1))
      m))

  (define (matrix-contract m1 m2 i k)
    (do ([j 0 (+ 1 j)]
         [accum 0 (+ accum (* (matrix-ref m1 i j)
                              (matrix-ref m2 j k)))])
        ((>= j (matrix-cols m1)) accum)))

  (define (matrix-map f m)
    (vector-map
     (lambda (v)
       (vector-map f v))
     m))

  ;; throw error if dims of a1 and a2 are bad
  ;; case on a1/a2 being a constant or not
  (define (mul a1 a2)
    (cond
     [(and (matrix? a1) (matrix? a2))
      (let* ([m (make-matrix (matrix-rows a1) (matrix-cols a2))])
        (do-matrix m m (i k) (matrix-set! m i k (matrix-contract a1 a2 i k))))]
     [(and (matrix? a2) (number? a1))
      (matrix-map (lambda (x) (* a1 x)) a2)]
     [(and (matrix? a1) (number? a2))
      (matrix-map (lambda (x) (* a2 x)) a1)]))

  (define (tr m)
    (do ([i 0 (+ 1 i)]
         [accum 0 (+ accum (matrix-ref m i i))])
        ((>= i (min (matrix-cols m) (matrix-rows m))) accum)))

  (define (T a1)
    (let ([m (make-matrix (matrix-cols a1) (matrix-rows a1))])
      (do-matrix m m (i j)
                 (matrix-set! m j i (matrix-ref a1 i j)))))

  ;; assert len(v) = len(w) and both vectors and all that
  (define (dot v w)
    (matrix-ref (mul (T v) w) 0))

  ;; vector->list probably inefficient, don't know how that's
  ;; implemented its probably O(n)
  (define (vector-fold-left f)
    (lambda (v)
      (fold-left f (vector-ref v 0) (vector->list v))))
  (define (matrix-fold-left f m)

    (fold-left f (matrix-ref m 0 0)
               (vector->list (vector-map (vector-fold-left f) m))))

  (define (matrix-fold f m)
    (let ([ret (matrix-ref m 0 0)])
      (do-matrix m ret (i j) (set! ret (f ret (matrix-ref m i j))))))

  ;; not good enough at macros to have matrix-map take variable args
  (define (matrix= A B)
    (matrix-fold-left
     (lambda (x y) (and x y))
     (vector-map
      (lambda (a b)
        (vector-map = a b))
      A B)))

  ;; Norms
  (define (matrix-max m)
    (matrix-fold-left max m))
  (define (matrix-min m)
    (matrix-fold-left min m))

  (define (frobenius-norm m)
    (matrix-fold-left
     +
     (matrix-map (lambda (x) (expt x 2)) m)))

  (define (euclidean-norm M)
    (cond
     [(matrix? M)
      (let ([ret 0])
        (matrix-map (lambda (x) (set! ret (+ ret (expt x 2)))) M))]
     [(vector? M)
      (do ([i 0 (+ 1 i)]
           [accum 0 (+ accum (expt (vector-ref M i) 2))])
          ((>= i (vector-length M)) (sqrt accum)))]))

  ;; Solve linear equation
  (define (linsolve input-A input-B)
    (when (and (matrix? input-A) (matrix? input-B))
      (let* ([A (matrix-copy input-A)]
             [B (matrix-copy input-B)]
             [nrows (matrix-rows A)]
             [ncols (matrix-cols A)])
        (let ([leading
               (lambda (col-idx)
                 (call/cc (lambda (break)
                            (do ([i col-idx (+ 1 i)])
                                ((or (>= col-idx nrows) (>= i ncols)) #f)
                              (when (not (zero? (matrix-ref A i col-idx)))
                                (break i))))))]
              [normalize-row
               (lambda (row col)
                 (let ([norm-val (matrix-ref A row col)])
                   (matrix-set-row! A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row A row)))
                   (matrix-set! B row 0
                                (/ (matrix-ref B row 0) norm-val))))]
              [swap-rows
               (lambda (row-idx1 row-idx2)
                 (unless (= row-idx1 row-idx2)
                   (let ([tmp-row (matrix-ref-row A row-idx1)]
                         [tmp-b (matrix-ref B row-idx1 0)])
                     (matrix-set-row! A row-idx1 (matrix-ref-row A row-idx2))
                     (matrix-set-row! A row-idx2 tmp-row)
                     (matrix-set! B row-idx1 0 (matrix-ref B row-idx2 0))
                     (matrix-set! B row-idx2 0 tmp-b))))]
              [eliminate-row
               (lambda (row-idx col-idx basis-idx)
                 (let ([ratio (matrix-ref A row-idx col-idx)])
                   (do ([i col-idx (+ 1 i)])
                       ((>= i ncols))
                     (matrix-set! A row-idx i
                                  (- (matrix-ref A row-idx i)
                                     (* ratio (matrix-ref A basis-idx i)))))
                   (matrix-set! B row-idx 0
                                (- (matrix-ref B row-idx 0)
                                   (* ratio (matrix-ref B basis-idx 0))))))])
          (call/cc
           (lambda (break)
             (do ([j 0 (+ 1 j)])
                 ((>= j ncols) B)
               (let ([lead (leading j)])
                 (when (not lead)
                   (break #f))
                 (normalize-row lead j)
                 (do ([i 0 (+ 1 i)])
                     ((>= i nrows))
                   (unless (= i lead)
                     (eliminate-row i j lead)))
                 (swap-rows j lead)))))))))

  (define (matrix-inverse input-A)
    (when (and (= (matrix-cols input-A) (matrix-rows input-A)) (matrix? input-A))
      (let* ([A (matrix-copy input-A)]
             [nrows (matrix-rows A)]
             [ncols (matrix-cols A)]
             [inv-A (matrix-identity nrows)])
        (let ([leading
               (lambda (col-idx)
                 (call/cc
                  (lambda (break)
                    (do ([i col-idx (+ 1 i)])
                        ((or (>= col-idx nrows) (>= i ncols)) #f)
                      (when (not (zero? (matrix-ref A i col-idx)))
                        (break i))))))]
              [normalize-row
               (lambda (row col)
                 (let ([norm-val (matrix-ref A row col)])
                   (matrix-set-row! A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row A row)))
                   (matrix-set-row! inv-A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row inv-A row)))))]
              [swap-rows
               (lambda (row-idx1 row-idx2)
                 (unless (= row-idx1 row-idx2)
                   (let ([tmp-row (matrix-ref-row A row-idx1)]
                         [tmp-row-inv (matrix-ref-row inv-A row-idx1)])
                     (matrix-set-row! A row-idx1 (matrix-ref-row A row-idx2))
                     (matrix-set-row! A row-idx2 tmp-row)
                     (matrix-set-row! inv-A row-idx1 (matrix-ref-row inv-A row-idx2))
                     (matrix-set-row! inv-A row-idx2 tmp-row-inv))))]
              [eliminate-row
               (lambda (row-idx col-idx basis-idx)
                 (let ([ratio (matrix-ref A row-idx col-idx)])
                   (do ([i 0 (+ 1 i)])
                       ((>= i ncols))
                     (matrix-set! A row-idx i
                                  (- (matrix-ref A row-idx i)
                                     (* ratio (matrix-ref A basis-idx i))))
                     (matrix-set! inv-A row-idx i
                                  (- (matrix-ref inv-A row-idx i)
                                     (* ratio (matrix-ref inv-A basis-idx i)))))))])
          (call/cc
           (lambda (break)
             (do ([j 0 (+ 1 j)])
                 ((>= j ncols) inv-A)
               (let ([lead (leading j)])
                 (when (not lead)
                   (break #f))
                 (normalize-row lead j)
                 (do ([i 0 (+ 1 i)])
                     ((>= i nrows))
                   (unless (= i lead)
                     (eliminate-row i j lead)))
                 (swap-rows j lead)))
             ))))))

  ;; Using the physics notation of the azimuthal 2\pi angle being phi
  (define (phi v)
    (atan (vector-ref v 1)
          (vector-ref v 0)))

  (define (theta v)
    (acos (/ (vector-ref v 2)
             (euclidean-norm v))))

  ;; assert dim 3
  (define (cross-product v w)
    (let ([cross-vw (make-vector 3)])
      (apply
       (lambda (ax ay az bx by bz)
         (vector-set! cross-vw 0 (- (* ay bz)
                                    (* by az)))
         (vector-set! cross-vw 1 (- (- (* ax bz)
                                       (* bx az))))
         (vector-set! cross-vw 2 (- (* ax by)
                                    (* bx ay))))
       (append (vector->list v) (vector->list w)))
      cross-vw))
  )
