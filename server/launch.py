"""Launch a local stdio MCP with project-local dependencies and write locations."""
import os
from pathlib import Path
import runpy
import site
import sys

sys.dont_write_bytecode=True
PACKAGE=Path(__file__).resolve().parent
WORKSPACE=Path(os.environ.get('MCP_WORKSPACE_ROOT',str(PACKAGE.parent))).resolve()
if WORKSPACE.drive.lower()=='c:':raise SystemExit('MCP workspace cannot be on C:')
if not PACKAGE.is_relative_to(WORKSPACE):raise SystemExit('Place the MCP package inside its configured project workspace')
os.environ['MCP_WORKSPACE_ROOT']=str(WORKSPACE)
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.environ['PYTHONIOENCODING']='utf-8'
temp=WORKSPACE/'docs/wiggle_mcp_runtime_tmp'
temp.mkdir(parents=True,exist_ok=True)
os.environ['TEMP']=str(temp);os.environ['TMP']=str(temp);os.environ['TMPDIR']=str(temp)
sys.path.insert(0,str(PACKAGE/'vendor'))
site.addsitedir(str(PACKAGE/'vendor'))
sys.path.insert(0,str(PACKAGE))
runpy.run_path(str(PACKAGE/'server.py'),run_name='__main__')
