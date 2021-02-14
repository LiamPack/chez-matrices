(library (matrix)
  (export make-matrix matrix-ref matrix-set! matrix-dims matrix-num-dims
          matrix-num-vals matrix-rows matrix-cols matrix? matrix-copy
          matrix-flatten vec->matrix matrix-ref-row matrix-set-row!
          matrix-contract matrix-ref-col matrix-set-col! matrix-map
          mul T matrix-fold-left tr euclidean-norm cross-product
          theta phi linsolve)
  (import (rnrs))

  (define (make-tensor dims init-val)
    (if (= (length dims) 1)
        (make-vector (car dims) init-val)
        (let ([rows (car dims)]
              [cols (cadr dims)])
          (do ([m (make-vector rows)]
               [i 0 (+ i 1)])
              ((= i rows) m)
            (vector-set! m i (make-tensor (cdr dims) init-val))))))

  (define (tensor-ref tensor indices)
    (let tensor-ref-rec ([tensor tensor] [indices (reverse indices)])
      (if (= (length indices) 1)
          (vector-ref tensor (car indices))
          (vector-ref (tensor-ref-rec tensor (cdr indices)) (car indices)))))

  (define (tensor-set! tensor indices val)
    (let ([base-vec (tensor-ref tensor (reverse (cdr (reverse indices))))])
      (vector-set! base-vec (car (reverse indices)) val)))

  (define (tensor-dims tensor)
    (let tensor-dims-rec ([tensor tensor] [dim-list (list)])
      (if (number? (vector-ref tensor 0))
          (reverse (cons (vector-length tensor) dim-list))
          (tensor-dims-rec (vector-ref tensor 0) (cons (vector-length tensor) dim-list)))))

  (define (tensor-rank tensor)
    (length (tensor-dims tensor)))

  (define-syntax do-nested
    (syntax-rules ()
      [(_ () e ...)  (begin e ...)]
      [(_ ((sym start end incr) i ...) e ...)
       (do ([sym start incr]) ((>= sym end) #f)
         (do-nested (i ...) e ...))]
      [(_ ((sym start end) i ...) e ...)
       (do ([sym start (+ 1 sym)]) ((>= sym end) #f)
         (do-nested (i ...) e ...))]
      [(_ ((sym end) i ...) e ...)
       (do ([sym 0 (+ 1 sym)]) ((>= sym end) #f)
         (do-nested (i ...) e ...))]))

  (define (make-matrix i j) (make-tensor (list i j) 0))
  (define (matrix-ref m i j) (tensor-ref m (list i j)))
  (define (matrix-set! m i j val) (tensor-set! m (list i j) val))
  (define matrix-dims tensor-dims)
  (define matrix-num-dims tensor-rank)
  (define (matrix-num-vals m) (map * (matrix-dims m)))
  (define (matrix-rows a) (vector-length a))
  (define (matrix-cols a)
    (if (vector? (vector-ref a 0))
        (vector-length (vector-ref a 0))
        1))

  (define matrix?
    (lambda (x)
      (and (vector? x)
           (> (vector-length x) 0)
           (vector? (vector-ref x 0)))))

  ;; assert proper dims
  (define (matrix-copy a)
    (let* ([m (make-matrix (matrix-rows a) (matrix-cols a))])
      (do-nested ([i (matrix-rows a)] [j (matrix-cols a)])
                 (matrix-set! m i j (matrix-ref a i j)))
      m))

  (define (matrix-copy! a b)
    (do-nested ([i (matrix-rows a)] [j (matrix-cols a)])
               (matrix-set! a i j (matrix-ref b i j))))

  (define (vec->matrix v)
    (make-matrix (vector-length v) 1))

  (define (matrix-flatten m)
    (let* ([v (make-vector (* (matrix-cols m) (matrix-rows m)))])
      (do-nested ([i (matrix-rows m)] [j (matrix-cols m)])
                 (vector-set! v (+ j (* i matrix-cols))
                              (matrix-ref m i j)))
      v))

  (define (matrix-ref-row m i) (vector-ref m i))
  (define (matrix-set-row! m i v) (vector-set! m i v))

  ;; verry inefficient but clean!
  (define (matrix-ref-col m i)
    (matrix-ref-row (T m) i))
  (define (matrix-set-col! m i v)
    (T (matrix-set-row! (T m) i v)))

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
     [(matrix? a2)
      (let* ([m (make-matrix (matrix-rows a1) (matrix-cols a2))])
        (do-nested ([i (matrix-rows a1)] [k (matrix-cols a2)])
                   (matrix-set! m i k (matrix-contract a1 a2 i k)))
        m)]
     [(number? a2)
      (matrix-map (lambda (x) (* a2 x)) a1)]))

  (define (tr m)
    (do ([i 0 (+ 1 i)]
         [accum 0 (+ accum (matrix-ref m i i))])
        ((>= i (min (matrix-cols m) (matrix-rows m))) accum)))

  (define (T a1)
    (let* ([m (make-matrix (matrix-cols a1) (matrix-rows a1))])
      (do-nested ([i 0 (matrix-rows a1)] [j 0 (matrix-cols a1)])
                 (matrix-set! m j i (matrix-ref a1 i j)))
      m))

  ;; assert len(v) = len(w) and both vectors and all that
  (define (dot v w)
    (matrix-ref (mul (T v) w)))

  ;; vector->list probably inefficient, don't know how that's
  ;; implemented its probably O(n)
  (define (vector-fold-left f)
    (lambda (v)
      (fold-left f (vector-ref v 0) (vector->list v))))
  (define (matrix-fold-left m f)
    (fold-left f (matrix-ref m 0 0)
               (vector->list (vector-map (vector-fold-left f) m))))

  ;; Norms
  (define (matrix-max m)
    (matrix-fold-left m max))
  (define (matrix-min m)
    (matrix-fold-left m min))

  (define (frobenius-norm m)
    (matrix-fold-left
     (matrix-map (lambda (x) (expt x 2)) m)
     +))

  (define (euclidean-norm M)
    (cond
     [(matrix? M) (frobenius-norm)]
     [(vector? M)
      (do ([i 0 (+ 1 i)]
           [accum 0 (+ accum (expt (vector-ref M i) 2))])
          ((>= i (vector-length M)) (sqrt accum)))]))

  ;; (define (frobebnius-norm)
  ;;   (sqrt (tr (mul m (T m)))))

  ;; SVD

  ;; Solve linear equation
  (define (linsolve input-A input-B)
    (when (and (matrix? input-A) (matrix? input-B))
      (let* ([A (matrix-copy input-A)]
             [B (matrix-copy input-B)]
             [nrows (matrix-rows A)]
             [ncols (matrix-cols A)])
        (letrec* ([leading
                   (lambda (col-idx)
                     (do ([i col-idx (+ 1 i)]
                          [ret #f])
                         ((or (>= col-idx nrows) (>= i ncols)) ret)
                       (when (not (zero? (matrix-ref A i col-idx)))
                         (set! ret i)
                         (set! i ncols))))]
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
                       (do ([i 0 (+ 1 i)])
                           ((>= i ncols) #f)
                         (let ([tmp (matrix-ref A row-idx1 i)])
                           (matrix-set! A row-idx1 i (matrix-ref A row-idx2 i))
                           (matrix-set! A row-idx2 i tmp)))
                       (let ([tmp (matrix-ref B row-idx1 0)])
                         (matrix-set! B row-idx1 0 (matrix-ref B row-idx2 0))
                         (matrix-set! B row-idx2 0 tmp))))]
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

          (do ([j 0 (+ 1 j)])
              ((>= j ncols))
            (let ([lead (leading j)])
              (if (not lead)
                  (set! j (+ 1 ncols)) ;;if i knew how to use
                  ;;call/cc it would happen
                  ;;here
                  (begin
                    (normalize-row lead j)
                    (do ([i 0 (+ 1 i)])
                        ((>= i nrows))
                      (unless (= i lead)
                        (eliminate-row i j lead)))
                    (swap-rows j lead)
                    ))))
          B))))

  ;; 3D vector operations

  ;; Assert vector-length(v) = 3

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
