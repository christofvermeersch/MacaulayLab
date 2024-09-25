% Create the document:
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{manual}[Software Manual Template]
\LoadClass[10pt,twoside]{article}

% Load all neccesary package:
\usepackage[paperheight = 24cm, paperwidth = 16cm,margin = 2cm]{geometry}
\usepackage{titlesec}
\usepackage{titleps}
\usepackage{titletoc}
\usepackage[english]{babel}
\usepackage{csquotes}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage{thmtools}
\usepackage{framed}
\usepackage{nicematrix}
\usepackage{pgfplotstable}
\usepackage{basictools}
\usepackage{multirow}
\usepackage{subcaption}
\usepackage[counterclockwise]{rotating}

% Initialize the Tikz/PgF environment:
\pgfplotsset{compat = 1.17, compat/bar nodes = 1.8}
\usetikzlibrary{patterns,decorations,calc,intersections,backgrounds}
\usepgfplotslibrary{fillbetween}


% Define the colors:



\definecolor{blue}{RGB}{69,117,180}
\definecolor{red}{RGB}{215,48,39}
\definecolor{yellow}{RGB}{255,199,100}
\definecolor{green}{RGB}{102,189,99}
\definecolor{purple}{RGB}{126,47,142}
\definecolor{string}{RGB}{181,49,247}
\definecolor{comment}{RGB}{0,128,19}

\colorlet{bluefill}{blue!25!white}
\colorlet{redfill}{red!25!white}
\colorlet{yellowfill}{yellow!25!white}
\colorlet{greenfill}{green!25!white}
\colorlet{purplefill}{purple!25!white}
\colorlet{grayfill}{gray!25!white}



% Create textual shorthands:
\newcommand{\maclab}{\textsf{MacaulayLab}}
\newcommand{\matlab}{\textsf{Matlab}}
\setlength{\fboxsep}{1pt}
\newcommand{\function}[1]{\colorbox{greenfill}{\texttt{\vphantom{phantom}#1}}}
\newcommand{\system}[1]{\texttt{#1}}
\newcommand{\problem}[1]{\texttt{#1}}
\newcommand{\toolbox}[1]{\textsf{#1}}

% Change cleveref options (to remove Eq. in front of equation number):
\crefname{equation}{}{}
\crefname{figure}{Figure}{Figures}

% Set the different units (lengths):
\newlength{\defaultthickness}
\setlength{\defaultthickness}{1pt}
\newlength{\ph}
\setlength{\ph}{\paperheight}
\newlength{\pw}
\setlength{\pw}{\paperwidth}
\newlength{\lw}
\setlength{\lw}{\linewidth}
\newlength{\la}
\setlength{\la}{4.7mm}
\newlength{\bh}
\setlength{\bh}{0.125\lw}
\newlength{\bw}
\setlength{\bw}{0.33\lw}
\newlength{\bs}
\setlength{\bs}{0.05\lw}
\newlength{\os}
\setlength{\os}{0.03\lw}
\newlength{\cu}
\setlength{\cu}{0.10\lw}
\newlength{\cwidth}
\setlength{\cwidth}{0.2\lw}
\newlength{\cheight}
\setlength{\cheight}{0.07\lw}
\newlength{\cskip}
\setlength{\cskip}{0.02\lw}
\setlength{\heavyrulewidth}{\defaultthickness}

% Change lay-out of the algorithm environment:
\patchcmd{\fs@ruled}
    {\def\@fs@post{\kern2pt\hrule\relax}}
    {\def\@fs@post{\kern\aboverulesep\hrule height\heavyrulewidth\relax}}
    {}{}
\floatstyle{ruled}
\restylefloat{algorithm}

% Define the different (text) environments:
\let\proof\relax

\newcommand{\codeinput}{>>\phantom{-}}
\newcommand{\codeoutput}{\\ \phantom{-} \\}
\newcommand{\codetilde}{\textasciitilde}
\newcommand{\tab}{\phantom{for-}}
\newcommand{\tabdash}{\phantom{for}}
\newcommand{\tabfor}{\phantom{>> }}

\newcommand{\createframe}[1]{
    \def\FrameCommand{\textcolor{#1}{\vrule width \heavyrulewidth} \hspace{5pt}}
    \MakeFramed {\advance\hsize - \width \FrameRestore}
}

\newenvironment{coderule}{\createframe{green}}{\endMakeFramed}
\declaretheoremstyle[preheadhook = \begin{coderule}, headfont = {\bfseries\color{green}}, notefont = {\bfseries\color{green}}, postfoothook = \end{coderule}\ignorespacesafterend\noindent]{codestyle}

\newenvironment{codeblock}{\colorlet{shadecolor}{greenfill} \ttfamily \begin{shaded*} \raggedright \noindent\ignorespaces}{\end{shaded*}}

\declaretheorem[name = Code, style = codestyle]{code}