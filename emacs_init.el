;; (c) Copyright Rosetta Commons Member Institutions.
;; (c) This file is part of the Rosetta software suite and is made available under license.
;; (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
;; (c) For more information, see http://www.rosettacommons.org. Questions about this can be
;; (c) addressed to University of Washington CoMotion, email: license@uw.edu.

;; @file rosetta_tests/emacs_init.el
;; @author Matthew O'Meara (mattjomeara@gmail.com)


;; This is a template emacs initialization file
;;
;; To use, copy all or some to  ~/.emacs.d/init.el


; Enable syntax hilighting
(global-font-lock-mode 1)


(setq inhibit-splash-screen t)

; short cut keys for compiling
(global-set-key [(f5)] 'compile)
(global-set-key [(f6)] 'recompile)

; like tab expansion for variable names etc.
(global-set-key [f12]         'dabbrev-expand)
(define-key esc-map [f12]     'dabbrev-completion)

; when two buffers are opened with the same name, they have the first
; unique subdirectory appeneded to the end of the buffer name to help
; distinguish the files
(require 'uniquify)
(setq uniquify-buffer-name-style 'reverse)
(setq uniquify-separator "/")
(setq uniquify-after-kill-buffer-p t) ; rename after killing uniquified
(setq uniquify-ignore-buffers-re "^\\*") ; don't muck with special buffers

; ido mode helps quickly open files, switch buffers etc.
; remember:
;   C-f -> go back to normal find-file
;   C-d -> open dired
;   C-j -> create new file (becuase enter would open the closest matching file)
(require 'ido)
(ido-mode t)
(setq ido-enable-flex-matching t)

; This function is broken on large python files like options_rosetta.py or pilot_apps.src.settings.all - it's a known bug https://github.com/bbatsov/prelude/issues/703).  Enable per-file with M-x which-function-mode, or turn this on here at your own risk and remember to disable it for large python files.
; put the current function in the bottom buffer
;(which-function-mode t)

; C-c o takes you from the cc file to hh file and visa versa
(add-hook 'c-mode-common-hook
  (lambda()
    (local-set-key  (kbd "C-c o") 'ff-find-other-file)))

;; Your init file should contain only one such instance.
;; If there is more than one, they won't work right.
(custom-set-variables
 '(safe-local-variable-values (quote ((rm-trailing-spaces . t) (show-trailing-whitespace . t) (rm-trailing-spaces . t)))))

; have tabbing in c++ mode work like it should with Rosetta
; if you find a place where it isn't working right
; do C-c C-s and it will tell you what syntax element you are in
; then add a c-set-offset command   -> '+ is one tab stop
(defun rosetta-c++-mode-hook ()
  ;; RosettaCommons Syntax style for C++ code

  (c-set-offset 'arglist-intro '+)
  (c-set-offset 'arglist-close 0)
  (c-set-offset 'arglist-cont-nonempty '+)
  (c-set-offset 'innamespace 0)
  (c-set-offset 'statement-cont '+)
  (c-set-offset 'stream-op '0)

 )
(add-hook 'c-mode-common-hook 'rosetta-c++-mode-hook)

; This hook deletes trailing whitespace, for any file type, in any mode
; upon saving.
; Uncomment to use.
;; (add-hook 'before-save-hook 'delete-trailing-whitespace)

; This function, and the below key redefinition, enable beautification of the current
; buffer by pressing the "F7" key
(defun rosetta-beautify (buffer-to-beautify)
  "Run the Rosetta beautification python script on the current file"
  (interactive)
  (if (eq (buffer-local-value 'major-mode buffer-to-beautify) 'c++-mode)
      ;Check to see if current buffer needs to be saved
      (progn
	(save-some-buffers nil (lambda () (eq (get-buffer (buffer-name)) buffer-to-beautify)))
	(if (not (buffer-modified-p buffer-to-beautify))
	    (progn
	      (with-current-buffer buffer-to-beautify (message (shell-command-to-string (concat "python ~/Rosetta/tools/python_cc_reader/beautifier.py --overwrite --filename " buffer-file-name))))
	      (with-current-buffer buffer-to-beautify (revert-buffer t t t))
	    )
	    (message "%s" (propertize "Not beautifying: buffer not saved" 'face '(:foreground "blue")))
	)
      )
      (message "%s" (propertize "Not beautifying: not currently in a c++-mode buffer" 'face '(:foreground "blue")))
      )
)
(global-set-key [(f7)] (lambda () (interactive) (rosetta-beautify (current-buffer))))
