set nocompatible
source $VIMRUNTIME/vimrc_example.vim
source $VIMRUNTIME/mswin.vim
behave mswin

filetype off

execute pathogen#infect()
execute pathogen#helptags()

filetype plugin indent on
syntax on

nnoremap <F4> :TagbarToggle<CR>
nnoremap <F8> :!ctags -R --python-kinds=-i *.py<CR>
nnoremap <F7> :NERDTreeToggle<CR>
nnoremap <F9> :%ArrangeColumn<CR>
nnoremap <F11> :%UnArrangeColumn<CR>


set gfn=Consolas:h10.5
set tabstop=4
set shiftwidth=4
set expandtab
set softtabstop=4
set smartindent
set autoindent
 
set hlsearch
set incsearch
set showmatch
 
set number
colorscheme desert
set background=dark

cd C:\Users\TZhu1\Documents

fu! SaveSess()
    execute 'mksession! ~/session.vim'
endfunction

fu! RestoreSess()
execute 'so ~/session.vim'
if bufexists(1)
    for l in range(1, bufnr('$'))
        if bufwinnr(l) == -1
            exec 'sbuffer ' . l
        endif
    endfor
endif
endfunction
" http://stackoverflow.com/questions/5142099/how-to-auto-save-vim-session-on-quit-and-auto-reload-on-start-including-split-wi

autocmd VimLeave * call SaveSess()
" autocmd VimEnter * call RestoreSess()
