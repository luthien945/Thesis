%
% University of Science and Technology of China Thesis Template
%
%   $Id: USTCthesis.cls 24 2011-05-16 09:29:01Z guolicai@gmail.com $
%
\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{USTCthesis}[2011/01/01 USTCthesis document class]

\newif\ifustc@master\ustc@masterfalse
\newif\ifustc@doctor\ustc@doctorfalse
\newif\ifustc@adobefont\ustc@adobefontfalse
\DeclareOption{master}{\ustc@mastertrue}
\DeclareOption{doctor}{\ustc@doctortrue}
\DeclareOption{adobefont}{\ustc@adobefonttrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{ctexbook}}
\ProcessOptions
\LoadClass[cs4size,a4paper,fancyhdr,fntef]{ctexbook}


\ifustc@master\relax\else
\ifustc@doctor\relax\else
\ClassError{USTCthesis}%
{You have to specify one of thesis options: master or doctor.}{}
\fi
\fi

\RequirePackage{times}
\RequirePackage{ccaption}
\RequirePackage[small,bf,labelsep=quad]{caption}
\RequirePackage{amsmath,amsthm,amsfonts,amssymb,bm}
\RequirePackage[lined,boxed,linesnumbered,algochapter]{algorithm2e}
\RequirePackage{graphicx,psfrag}
\RequirePackage{titlesec}
\RequirePackage{subfig}
\RequirePackage{titletoc}
\RequirePackage[left=3.17cm,right=3.17cm,top=2.54cm,bottom=2.54cm]{geometry}%left=3.17cm,right=3.17cm,top=2.54cm,bottom=2.54cm,,headheight=0.5cm,
\RequirePackage{enumerate}
\RequirePackage{multirow}
\RequirePackage{subfig,epsfig}
\RequirePackage{float}
\RequirePackage{paralist}
\RequirePackage{booktabs}
\RequirePackage{setspace}
\RequirePackage[perpage,symbol]{footmisc}
\RequirePackage{pifont}
\RequirePackage{hyperref}
\RequirePackage[absolute]{textpos}
\RequirePackage[numbers,sort&compress,square]{natbib}

    
\hypersetup
{
    CJKbookmarks=true,
    bookmarksnumbered=true,
    bookmarksopen=true,
    breaklinks=true,
    colorlinks=true,
    linkcolor=black,
    citecolor=black,
    plainpages=false,
    pdfstartview=FitH
}

\ifustc@adobefont
\setCJKmainfont[BoldFont={Adobe Heiti Std}, ItalicFont={Adobe Kaiti Std}]{Adobe Song Std}
\setCJKsansfont{Adobe Heiti Std}
\setCJKmonofont{Adobe Fangsong Std}

\setCJKfamilyfont{song}{Adobe Song Std}
\setCJKfamilyfont{hei}{Adobe Heiti Std}
\setCJKfamilyfont{kai}{Adobe Kaiti Std}
\setCJKfamilyfont{fs}{Adobe Fangsong Std}
\setCJKfamilyfont{verthei}[RawFeature={vertical:}]{Adobe Heiti Std}%{Adobe Heiti Std} 这里还是改成Adobe Heiti Std吧，在Linux下面SimHei貌似不行，改成这个就好了
\else
\setCJKmainfont[BoldFont={SimHei}, ItalicFont={KaiTi}]{NSimSun}
\setCJKsansfont{SimHei}
\setCJKmonofont{FangSong}

\setCJKfamilyfont{song}{NSimSun}
\setCJKfamilyfont{hei}{SimHei}
\setCJKfamilyfont{kai}{KaiTi}   % XP对应 KaiTi_GB2312，Vista对应KaiTi，注意根据系统切换
\setCJKfamilyfont{fs}{FangSong} % XP对应 FangSong_GB2312，Vista对应FangSong，注意根据系统切换
\setCJKfamilyfont{verthei}[RawFeature={vertical:}]{SimHei}%{Adobe Heiti Std}
\fi

\newcommand{\song}{\CJKfamily{song}}    % 宋体
\newcommand{\fs}{\CJKfamily{fs}}        % 仿宋体
\newcommand{\kai}{\CJKfamily{kai}}      % 楷体
\newcommand{\hei}{\CJKfamily{hei}}      % 黑体
\newcommand{\verthei}{\CJKfamily{verthei}}

\setmainfont{Times New Roman}
\setsansfont[BoldFont={Arial:style=Bold}]{Arial}
\setmonofont[BoldFont={Courier New Bold}]{Courier New}

%XeTex中文断行
\XeTeXlinebreaklocale "zh"
\XeTeXlinebreakskip = 0pt plus 1pt minus 0.1pt

%titleformat
\setcounter{secnumdepth}{3}

\CTEXsetup[
    number={\arabic{chapter}},
    format={\centering},
    nameformat={\hei\fontsize{16pt}{16pt}\linespread{1}\selectfont\bfseries},
    titleformat={\hei\fontsize{16pt}{16pt}\linespread{1}\selectfont\bfseries},
    beforeskip={14pt},
    afterskip={42pt} % 与节的skip重叠时选较大的skip，章标题段后18pt+节标题24pt
]{chapter}

\CTEXsetup[
    format={\flushleft\hei\fontsize{14pt}{14pt}\linespread{1}\selectfont},
    aftername={\quad},
    beforeskip={24pt},
    afterskip={6pt}
]{section}

\CTEXsetup[
    format={\flushleft\hei\fontsize{13pt}{13pt}\linespread{1}\selectfont},
    aftername={\quad},
    beforeskip={12pt},
    afterskip={6pt}
]{subsection}

\CTEXsetup[
    format={\flushleft\hei\fontsize{12pt}{12pt}\linespread{1}\selectfont},
    aftername={\quad},
    beforeskip={10pt},
    afterskip={6pt}
]{subsubsection}

 \captionsetup[figure]
{
    aboveskip=6pt,
    belowskip=12pt,
}

\captionsetup[tabel]
{
    aboveskip=6pt,
    belowskip=6pt,
}


\newtheorem{theorem}{\textbf{\hspace{0.7cm}定理}}[section]
\newtheorem{lemma}{\textbf{\hspace{0.7cm}引理}}[section]
\newtheorem{example}{\textbf{\hspace{0.7pt}例}}[section]
\newtheorem{definition}{\textbf{\hspace{0.7pt}定义}}[section]
\newtheorem{axiom}{\textbf{\hspace{0.7pt}公理}}[section]
\newtheorem{property}{\textbf{\hspace{0.7pt}性质}}[section]
\newtheorem{proposition}{\textbf{\hspace{0.7pt}命题}}[section]
\newtheorem{corollary}{\textbf{\hspace{0.7pt}推论}}[section]
\newtheorem{remark}{\textbf{\hspace{0.7pt}注解}}[section]
\newtheorem{condition}{\textbf{\hspace{0.7pt}条件}}[section]
\newtheorem{conclusion}{\textbf{\hspace{0.7pt}结论}}[section]
\newtheorem{assumption}{\textbf{\hspace{0.7pt}假设}}[section]
\newtheorem{prove}{\textbf{\hspace{0.7cm}证明}}[section]
\renewcommand{\listalgorithmcfname}{算法}
\renewcommand{\algorithmcfname}{算法}
%\renewcommand\appendixname{附录}

%行距，段落设置
\renewcommand{\normalsize}
    {\fontsize{12pt}{20pt}\linespread{1}\selectfont}   %默认字体12pt，行距20pt
%\renewcommand{\small}
%    {\fontsize{10.5pt}{15.6pt}\linespread{1}\selectfont}
    %小号字体10.5pt，行距15.6pt(五号的单倍行距)

\setlength{\parskip}{0pt plus3pt minus0pt}%{1ex plus1.5ex minus0.5ex}%{0.5\baselineskip}%{1.5pt plus1pt minus1pt}

\newcommand\ustc@underline[2][6em]{\hskip1pt\underline{\hb@xt@ #1{\hss#2\hss}}\hskip3pt}
\def\cleardoublepage
{
    \clearpage
    \if@twoside 
        \ifodd\c@page
        \else
            \thispagestyle{empty}%
            \hbox{}\newpage
            \if@twocolumn\hbox{}\newpage
            \fi
        \fi
    \fi
}

\renewcommand{\floatpagefraction}{0.80}
\newcommand\USTCTspace{\protect\CTEX@spaceChar}

\def\ustc@define@term#1{
    \expandafter\gdef\csname #1\endcsname##1{%
        \expandafter\gdef\csname ustc@#1\endcsname{##1}}
    \csname #1\endcsname{}}

\ustc@define@term{depart}
\ustc@define@term{degree}
\renewcommand\title[1]{\def\ustc@title{#1}}
\renewcommand\author[1]{\def\ustc@author{#1}}
\ustc@define@term{advisor}
\ustc@define@term{coadvisor}
\ustc@define@term{major}
\ustc@define@term{submitdate}

\renewcommand\maketitle{%
    \pdfbookmark[-1]{\ustc@title}{title}
    \cleardoublepage
    \thispagestyle{empty}
    \setlength{\TPHorizModule}{1mm}
    \setlength{\TPVertModule}{1mm}

    \begin{textblock}{150}[0.5,0](105,49)
        \begin{center}
            \includegraphics[scale=1]{ustc_logo_text.eps}
        \end{center}
    \end{textblock}

     \begin{textblock}{150}[0.5,0.5](105,85)
        \begin{center}
        \ifustc@doctor
            \hei\fontsize{56}{67.2}\selectfont {博士学位论文}
        \else
            \hei\fontsize{56}{67.2}\selectfont {硕士学位论文}
         \fi
        \end{center}
    \end{textblock}
   
    \begin{textblock}{150}[0.5,0.5](105,135)
        \begin{center}
        \includegraphics[width=4.52cm,height=4.52cm]{ustc_logo_fig.eps}
        \end{center}
    \end{textblock}
      
    \begin{textblock}{120}[0.5,0.613](105,183)
        \begin{center}
        \parbox[c]{11cm}
            {\center\hei\fontsize{26}{33}\selectfont\bfseries\ustc@title}
        \end{center}
    \end{textblock}

    \begin{textblock}{120}[0.5,1](105,260)
    \begin{center}
        { \fontsize{16}{21.5}\selectfont
        \begin{tabular}{p{3.8cm}p{7.5cm}}
            \hei{作者姓名：} & \hfil\ustc@author\hfil\\
            \hei{学科专业：} &\hfil\ustc@major\hfil\\
            \hei{导师姓名：} & \hfil\ustc@advisor\hfil\\
            \hei{完成时间：} &\hfil\ustc@submitdate\hfil 
        \end{tabular}}
    \end{center}
    \end{textblock}
    \null
    \clearpage
    \if@twoside
        \thispagestyle{empty}
        \cleardoublepage
    \fi
}

\ustc@define@term{etitle}
\ustc@define@term{eauthor}
\ustc@define@term{eadvisor}
\ustc@define@term{ecoadvisor}
\ustc@define@term{esubmitdate}
\ustc@define@term{emajor}

\newcommand\makeetitle{%
    %\pdfbookmark{英文标题}{etitle}
    \cleardoublepage
    \thispagestyle{empty}
    \setlength{\TPHorizModule}{1mm}
    \setlength{\TPVertModule}{1mm}

    \begin{textblock}{150}[0.5,0](105,47.5)
        \begin{center}
         \fontsize{20}{30}\selectfont\textsf{University of Science and Technology of China}
        \end{center}
    \end{textblock}

    \begin{textblock}{150}[0.5,0](105,56.5)
        \begin{center}
        \ifustc@doctor
            \fontsize{26}{30}\selectfont\textsf{A dissertation for doctor degree}
        \else
            \fontsize{26}{30}\selectfont\textsf{A dissertation for master degree}
        \fi
        \end{center}
    \end{textblock}

    \begin{textblock}{150}[0.5,0](105,88)
        \begin{center}
        \includegraphics[width=4.52cm,height=4.52cm]{ustc_logo_fig.eps}
        \end{center}
    \end{textblock}

    \begin{textblock}{110}[0.5,0.613](105,183)
        \begin{center}
            \fontsize{26}{30}\selectfont\bfseries\textsf{\ustc@etitle}
        \end{center}
    \end{textblock}

    \begin{textblock}{120}[0.5,1](104,260)
        \begin{center}
        { \fontsize{16pt}{21.5pt}\selectfont
        \begin{tabular}{p{3.8cm}p{7.5cm}}
            Author's Name :& \hfil\ustc@eauthor\hfil\\
            Speciality :& \hfil\ustc@emajor\hfil\\
            Supervisor :& \hfil\ustc@eadvisor\hfil\\
            Finished Time :& \hfil\ustc@esubmitdate\hfil
        \end{tabular}}
        \end{center}
    \end{textblock}
    \null
    \clearpage
    \if@twoside
        \thispagestyle{empty}
        \cleardoublepage
    \fi
}

\newcommand\makebookspine{
    %\pdfbookmark{书\hspace{1em}脊}{bookmark}
    \cleardoublepage
        \thispagestyle{empty}
    \begin{center}
    \rotatebox{-90}{
        \framebox[0.97\textheight][c]{\ziju{0.17}\verthei\zihao{-4} \raisebox{10pt}{ \ustc@title
            \hspace*{20pt}  \ustc@depart \hspace*{20pt}  \ustc@author  \hspace*{30pt}{中国科学技术大学} }}
    }
    \end{center}
    \clearpage
        \if@twoside
        \thispagestyle{empty}
    \cleardoublepage
        \fi
}

%% 定义页眉横线
\newcommand{\makeheadrule}{%
    %\rule[.7\baselineskip]{\headwidth}{1.2pt}\vskip-0.97\baselineskip%
     \rule[.4\baselineskip]{\headwidth}{0.75pt}}%

%\makeatletter
\renewcommand{\headrule}{%
    {\if@fancyplain\let\headrulewidth\plainheadrulewidth\fi}%
    \makeheadrule}
%\makeatother

\newcommand{\ustc@declaretext}{%
    本人声明所呈交的学位论文,是本人在导师指导下进行研究工作所取得的成果。
    除已特别加以标注和致谢的地方外，论文中不包含任何他人已经发表或撰写过
    的研究成果。与我一同工作的同志对本研究所做的贡献均已在论文中作了明确的说明。
}

\newcommand{\ustc@authorization}{%
    作为申请学位的条件之一，学位论文著作权拥有者授权中国科学技术大学拥有
    学位论文的部分使用权，即：学校有权按有关规定向国家有关部门或机构送交
    论文的复印件和电子版，允许论文被查阅和借阅，可以将学位论文编入《中国
    学位论文全文数据库》等有关数据库进行检索，可以采用影印、缩印或扫描等
    复制手段保存、汇编学位论文。本人提交的电子文档的内容和纸质论文的内容
    相一致。

    保密的学位论文在解密后也遵守此规定。
}

\def\ustc@secret@newtitle{□\ 公开\hspace{0.5cm}□\ 保密\ \ustc@underline[4em]\relax\hspace{0.1em}年}

\newcommand\makeauthorization{%
    % \pdfbookmark{原创与授权声明}{authorization}
    \cleardoublepage
    \thispagestyle{empty}
    \begin{center}\zihao{3}\hei{中国科学技术大学学位论文原创性声明}\end{center}
    \par\ustc@declaretext\par
    \vskip 1.0cm
    {\begin{tabular}{p{6cm}p{6cm}}
        作者签名：\ustc@underline{} & 签字日期：\ustc@underline{}
     \end{tabular}}
     \vskip2.0cm
     \begin{center}\zihao{3}\hei{中国科学技术大学学位论文授权使用声明}\end{center}
     \par\ustc@authorization\par
     \vskip0.5cm
     \ustc@secret@newtitle\par
     \vskip1.0cm
     {\renewcommand{\arraystretch}{2.0}%行高
             \begin{tabular}{p{6cm}p{6cm}}
                 作者签名：\ustc@underline{} & 导师签名：\ustc@underline{}\\
                 签字日期：\ustc@underline{} & 签字日期：\ustc@underline{}
             \end{tabular}}
     \clearpage
     \if@twoside
        \thispagestyle{empty}
        \cleardoublepage
     \fi
}

\fancypagestyle{plain}{}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[C]{\small \leftmark}

\renewcommand{\headrulewidth}{0.75pt}

\renewcommand\frontmatter{%
    \cleardoublepage
    \@mainmatterfalse
    \pagenumbering{Roman}
    \fancyfoot{}
    \fancyfoot[C]{\small \thepage}
}

\renewcommand\mainmatter{%
    \cleardoublepage
    \@mainmattertrue
    \pagenumbering{arabic}
    \fancyfoot{}
    \fancyfoot[RO,LE]{\small \thepage}
}

\newenvironment{abstract}
{\chapter{摘\hspace{1em}要}}
{}
\newcommand\keywords[1]{%
    \vspace{3.5ex}\noindent{\hei 关键词：} #1}

\newenvironment{eabstract}
{
    \chapter{\textsf{ABSTRACT}}
    \markboth{Abstract}{Abstract}
}
{ }

\newcommand\ekeywords[1]{%
    \vspace{3.5ex}\noindent{\bfseries {Keywords:~}} #1}

    %tableofcontents
    \renewcommand{\contentsname}{目\hspace{1em}录}
    \renewcommand\tableofcontents{%
        \if@twocolumn
            \@restonecoltrue\onecolumn
            \else
            \@restonecolfalse
            \fi
            \chapter{\contentsname}%
            \renewcommand{\baselinestretch}{1}
        \@mkboth{\MakeUppercase\contentsname}{\MakeUppercase\contentsname}%
            \@starttoc{toc}%
            \if@restonecol\twocolumn\fi
    }
%% tableofcontent
\titlecontents{chapter}[0em]
{\vspace{6bp}\hei\zihao{-4}}
{\thecontentslabel\hspace*{0.5em}}
{}
{\titlerule*[0.5pc]{$\cdot$}\contentspage}
\titlecontents{section}[1.1\ccwd]
{\vspace{1bp} \zihao{-4}}
{\thecontentslabel\hspace*{0.5em}}
{}
{\titlerule*[0.5pc]{$\cdot$}\contentspage}
\titlecontents{subsection}[2.3\ccwd]
{\vspace{1bp} \zihao{-4}}
{\thecontentslabel\hspace*{0.5em}}
{}
{\titlerule*[0.5pc]{$\cdot$}\contentspage}


\renewcommand{\bibsection}{
    \chapter*{\bibname}
    \markboth{参考文献}{参考文献}
    \addcontentsline{toc}{chapter}{参考文献}
    \fontsize{10.5pt}{20pt}\selectfont
    \@mkboth{\MakeUppercase\bibname}{\MakeUppercase\bibname}
    \setlength{\bibsep}{0ex}
}

\renewenvironment{thanks}
{
    \chapter*{致$\quad$谢}
    \addcontentsline{toc}{chapter}{致谢}  
    \markboth{致$\quad$谢}{致$\quad$谢}
}
{}

\newenvironment{denotation}[1][2.5cm]{
    \chapter{主要符号对照表} % no tocline
        \noindent\begin{list}{}%
        {\vskip-30bp%\xiaosihao[1.6]
            \renewcommand\makelabel[1]{##1\hfil}
            \setlength{\labelwidth}{#1} % 标签盒子宽度
            \setlength{\labelsep}{0.5cm} % 标签与列表文本距离
            \setlength{\itemindent}{0cm} % 标签缩进量
            \setlength{\leftmargin}{\labelwidth+\labelsep} % 左边界
            \setlength{\rightmargin}{0cm}
            \setlength{\parsep}{0cm} % 段落间距
            \setlength{\itemsep}{0cm} % 标签间距
            \setlength{\listparindent}{0cm} % 段落缩进量
            \setlength{\topsep}{0pt} % 标签与上文的间距
        }}
    {\end{list}
}
\endinput
%%
%% End of file `USTCthesis.cls'.
