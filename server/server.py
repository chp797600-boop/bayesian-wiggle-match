"""MCP stdio tools. The launcher loads the isolated official SDK first."""
import os
from pathlib import Path
from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from app import Backend

WORKSPACE=Path(os.environ['MCP_WORKSPACE_ROOT']).resolve()
backend=Backend(WORKSPACE)
mcp=FastMCP('Bayesian annual tree-ring wiggle matching',instructions=(
    'Use measured single-year reference segments with published independent calendar dates. '
    'List references and inspect provenance before selecting one. Explicit reference, search and offset prior required. '
    'Library gaps and years outside coverage are unassessed, not excluded. Do not pool references or claim automatic reliable dating. '
    'All calendar bounds are astronomical integer years: 0=1 BCE. Ages and offsets are conventional 14C BP.'))


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True))
def list_reference_data(start_astro_year:int|None=None,end_astro_year:int|None=None,include_demo:bool=False)->dict:
    """List bundled annual 14C reference segments, coverage gaps and calendar/source evidence. Date overlap alone is not suitability."""
    return backend.list_reference_data(start_astro_year,end_astro_year,include_demo)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True))
def get_reference_info(reference_id:str)->dict:
    """Verify a bundled reference hash and return provenance, physical-tree aliases and limitations."""
    return backend.get_reference_info(reference_id)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True))
def validate_floating_series(samples_path:str)->dict:
    """Check a project CSV with sample_id, ring_offset, c14_age_bp, sigma_bp, sample_span_years, measurement_ids, tree_id. Not scientific acceptance."""
    return backend.validate_floating_series(samples_path)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False,destructiveHint=False))
def run_wiggle_match(samples_path:str,reference_id:str,start_min:int,start_max:int,
                     offset_mean_bp:float,offset_sd_bp:float,run_id:str,allow_incomplete:bool=False)->dict:
    """Match to one explicitly selected measured reference; write a new project output folder. Prior SD is explicit; never choose it to fit truth. Missing coverage requires explicit exploratory opt-in."""
    return backend.run_wiggle_match(samples_path,reference_id,start_min,start_max,
                                    offset_mean_bp,offset_sd_bp,run_id,allow_incomplete)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False,destructiveHint=False))
def assess_match(run_id:str,assessment_id:str,seed:int,draws:int=256)->dict:
    """Run saved-input verification, leave-one-ring-out and posterior predictive checks. Does not supply final shared-offset Acomb or automatic five-gate acceptance."""
    return backend.assess_match(run_id,assessment_id,seed,draws)


if __name__=='__main__':mcp.run(transport='stdio')
