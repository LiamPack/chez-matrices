(library (matrix)
  (export make-array array-ref array-set! array-dims array-num-dims
          array-num-vals array-rows array-cols matrix? array-copy
          array-flatten vec->array array-ref-row array-set-row!
          array-contract array-ref-col array-set-col! array-map
          mul T array-fold-left tr)
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

  (define (make-array i j) (make-tensor (list i j) 0))
  (define (array-ref m i j) (tensor-ref m (list i j)))
  (define (array-set! m i j val) (tensor-set! m (list i j) val))
  (define array-dims tensor-dims)
  (define array-num-dims tensor-rank)
  (define (array-num-vals m) (map * (array-dims m)))
  (define (array-rows a) (vector-length a))
  (define (array-cols a)
    (if (vector? (vector-ref a 0))
        (vector-length (vector-ref a 0))
        1))

  (define matrix?
    (lambda (x)
      (and (vector? x)
           (> (vector-length x) 0)
           (vector? (vector-ref x 0)))))

  ;; assert proper dims
  (define (array-copy a)
    (let* ([m (make-array (array-cols a) (array-rows a))])
      (do-nested ([i (array-rows a)] [j (array-cols a)])
                 (array-set! m i j (array-ref a i j)))
      m))

  (define (vec->array v)
    (make-array (vector-length v) 1))

  (define (array-flatten m)
    (let* ([v (make-vector (* (array-cols m) (array-rows m)))])
      (do-nested ([i (array-rows m)] [j (array-cols m)])
                 (vector-set! v (+ j (* i array-cols))
                              (array-ref m i j)))
      v))

  (define (array-ref-row m i) (vector-ref m i))
  (define (array-set-row! m i v) (vector-set! m i v))

  ;; verry inefficient but clean!
  (define (array-ref-col m i)
    (array-ref-row (T m) i))
  (define (array-set-col! m i v)
    (T (array-set-row! (T m) i v)))

  (define (array-contract m1 m2 i k)
    (do ([j 0 (+ 1 j)]
         [accum 0 (+ accum (* (array-ref m1 i j)
                              (array-ref m2 j k)))])
        ((>= j (array-cols m1)) accum)))

  (define (array-map f m)
    (vector-map
     (lambda (v)
       (vector-map f v))
     m))

  ;; throw error if dims of a1 and a2 are bad
  ;; case on a1/a2 being a constant or not
  (define (mul a1 a2)
    (cond
     [(matrix? a2)
      (let* ([m (make-array (array-rows a1) (array-cols a2))])
        (do-nested ([i (array-rows a1)] [k (array-cols a2)])
                   (array-set! m i k (array-contract a1 a2 i k)))
        m)]
     [(number? a2)
      (array-map (lambda (x) (* a2 x)) a1)]))

  (define (tr m)
    (do ([i 0 (+ 1 i)]
         [accum 0 (+ accum (array-ref m i i))])
        ((>= i (min (array-cols m) (array-rows m))) accum)))

  (define (T a1)
    (let* ([m (make-array (array-cols a1) (array-rows a1))])
      (do-nested ([i 0 (array-rows a1)] [j 0 (array-cols a1)])
                 (array-set! m j i (array-ref a1 i j)))
      m))

  ;; assert len(v) = len(w) and both vectors and all that
  (define (dot v w)
    (array-ref (mul (T v) w)))

  ;; vector->list probably inefficient, don't know how that's
  ;; implemented its probably O(n)
  (define (vector-fold-left f)
    (lambda (v)
      (fold-left f (vector-ref v 0) (vector->list v))))
  (define (array-fold-left m f)
    (fold-left f (array-ref m 0 0)
               (vector->list (vector-map (vector-fold-left f) m))))

  ;; Norms
  (define (array-max m)
    (array-fold-left m max))
  (define (array-min m)
    (array-fold-left m min))

  (define (frobebnius-norm)
    (sqrt (tr (mul m (T m)))))

  ;; SVD

  ;; Solve linear equation
  )
