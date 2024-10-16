% Create the document:
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{manual}[Software Manual Template]
\LoadClass[10pt,twoside]{article}

% Load all neccesary package:
\usepackage[paperheight = 24cm, paperwidth = 16cm,margin = 2cm]{geometry}
\usepackage[english]{babel}
\usepackage{csquotes}
\usepackage{amsthm}
\usepackage{thmtools}
\usepackage{framed}
\usepackage{xcolor}
\usepackage{tikz}
\usepackage{pgfplots}
\usepackage{siunitx}
\usepackage{booktabs}
\usepackage{mathtools}
\usepackage{bm}
\usepackage{float}
\usepackage[colorlinks]{hyperref}
\usepackage{cleveref}

% Initialize the Tikz/PgF environment:
\pgfplotsset{compat = 1.17, compat/bar nodes = 1.8}
\usetikzlibrary{patterns,decorations,calc,intersections,backgrounds}

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
\newcommand{\function}[1]{\colorbox{greenfill}{\texttt{\vphantom{phantom}#1}}}
\newcommand{\system}[1]{\texttt{#1}}
\newcommand{\problem}[1]{\texttt{#1}}
\newcommand{\toolbox}[1]{\textsf{#1}}
\newcommand{\vc}[1]{\bm{#1}}
\newcommand{\mt}[1]{\bm{#1}}
\newcommand{\mcl}[1]{\bm{\mathcal{#1}}}

% Change cleveref options:
\crefname{equation}{}{}
\crefname{figure}{Figure}{Figures}

% Set the different units (lengths):
\newlength{\defaultthickness}
\setlength{\defaultthickness}{1pt}
\newlength{\cu}
\setlength{\cu}{0.10\linewidth}
\setlength{\heavyrulewidth}{\defaultthickness}
\setlength{\fboxsep}{1pt}

% Change lay-out of the algorithm environment:
\patchcmd{\fs@ruled}
    {\def\@fs@post{\kern2pt\hrule\relax}}
    {\def\@fs@post{\kern\aboverulesep\hrule height\heavyrulewidth\relax}}
    {}{}
\floatstyle{ruled}
\restylefloat{algorithm}

% Define the code environment:
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