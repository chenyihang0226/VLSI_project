param(
    [int]$Threads = 0,
    [int]$K = 43
)

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $ScriptDir

try {
    $Exe = Join-Path $ScriptDir "main.exe"
    $BenchDir = Join-Path $ScriptDir "ICCAD2021-TopoPart-Benchmarks\Generated Benchmarks"
    $TopoFile = Join-Path $ScriptDir "ICCAD2021-TopoPart-Benchmarks\FPGA Graph\MFS2"
    $ResultDir = Join-Path $ScriptDir "results"
    $LogDir = Join-Path $ResultDir "mfs2_logs"
    $SummaryFile = Join-Path $ResultDir "mfs2_summary.csv"

    if (!(Test-Path $Exe)) {
        throw "main.exe not found: $Exe"
    }
    if (!(Test-Path $BenchDir)) {
        throw "Benchmark directory not found: $BenchDir"
    }
    if (!(Test-Path $TopoFile)) {
        throw "Topology file not found: $TopoFile"
    }

    New-Item -ItemType Directory -Force -Path $ResultDir | Out-Null
    New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

    "case,exit_code,wall_seconds,cutsize,algo_time_seconds,log_file" |
        Set-Content -Path $SummaryFile -Encoding UTF8

    $Cases = Get-ChildItem -Path $BenchDir -File |
        Where-Object { $_.Name -match '^case\d+$' } |
        Sort-Object {
            if ($_.Name -match '^case(\d+)$') { [int]$Matches[1] } else { 0 }
        }

    Write-Host "========================================"
    Write-Host "MFS2 topology benchmark test"
    Write-Host "Working directory: $ScriptDir"
    Write-Host "Benchmark dir: $BenchDir"
    Write-Host "Topology: $TopoFile"
    Write-Host "Threads: $Threads"
    Write-Host "K: $K"
    Write-Host "Summary: $SummaryFile"
    Write-Host "========================================"

    foreach ($Case in $Cases) {
        $CaseName = $Case.Name
        $LogFile = Join-Path $LogDir "$CaseName.log"
        $Start = Get-Date

        Write-Host ""
        Write-Host "========================================"
        Write-Host "Running $CaseName"
        Write-Host "Command: $Exe `"$($Case.FullName)`" $Threads $K --topo `"$TopoFile`""
        Write-Host "Log: $LogFile"
        Write-Host "========================================"

        $StdoutFile = Join-Path $LogDir "$CaseName.stdout.tmp"
        $StderrFile = Join-Path $LogDir "$CaseName.stderr.tmp"
        Remove-Item -LiteralPath $StdoutFile, $StderrFile -Force -ErrorAction SilentlyContinue

        $Args = @(
            "`"$($Case.FullName)`"",
            "$Threads",
            "$K",
            "--topo",
            "`"$TopoFile`""
        )
        $Process = Start-Process -FilePath $Exe `
            -ArgumentList $Args `
            -Wait `
            -NoNewWindow `
            -PassThru `
            -RedirectStandardOutput $StdoutFile `
            -RedirectStandardError $StderrFile
        $ExitCode = $Process.ExitCode

        $Output = @()
        if (Test-Path $StdoutFile) {
            $Output += Get-Content -Path $StdoutFile
        }
        if (Test-Path $StderrFile) {
            $Output += Get-Content -Path $StderrFile
        }
        $Output | Tee-Object -FilePath $LogFile
        Remove-Item -LiteralPath $StdoutFile, $StderrFile -Force -ErrorAction SilentlyContinue

        $End = Get-Date
        $WallSeconds = [math]::Round(($End - $Start).TotalSeconds, 3)
        $Cutsize = ""
        $AlgoTime = ""

        foreach ($Line in $Output) {
            $Text = [string]$Line
            if ($Text -match 'Best cutsize found:\s*(\d+)') {
                $Cutsize = $Matches[1]
            }
            if ($Text -match 'Cutsize:\s*(\d+)\s+Time:\s*([0-9.]+)\s*s') {
                $Cutsize = $Matches[1]
                $AlgoTime = $Matches[2]
            }
        }

        "$CaseName,$ExitCode,$WallSeconds,$Cutsize,$AlgoTime,$LogFile" |
            Add-Content -Path $SummaryFile -Encoding UTF8

        Write-Host "Finished ${CaseName}: exit=$ExitCode, wall=${WallSeconds}s, cutsize=$Cutsize, algo_time=${AlgoTime}s"
    }

    Write-Host ""
    Write-Host "========================================"
    Write-Host "All MFS2 benchmark tests finished."
    Write-Host "Summary: $SummaryFile"
    Write-Host "Logs: $LogDir"
    Write-Host "========================================"
}
finally {
    Pop-Location
}
